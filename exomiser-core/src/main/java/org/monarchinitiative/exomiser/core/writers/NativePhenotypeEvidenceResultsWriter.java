package org.monarchinitiative.exomiser.core.writers;

import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVPrinter;
import org.monarchinitiative.exomiser.core.analysis.AnalysisResults;
import org.monarchinitiative.exomiser.core.prioritisers.HiPhivePriorityResult;
import org.monarchinitiative.exomiser.core.prioritisers.PhivePriorityResult;
import org.monarchinitiative.exomiser.core.prioritisers.model.GeneModelPhenotypeMatch;
import org.monarchinitiative.exomiser.core.prioritisers.model.GeneOrthologModel;
import org.monarchinitiative.exomiser.core.phenotype.TransdCheckpoint;
import org.monarchinitiative.exomiser.core.phenotype.ScoringResourceDigest;

import java.io.IOException;
import java.io.StringWriter;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.List;

/** Auditable raw model compatibility rows for native phenotype scorer substitution. */
public final class NativePhenotypeEvidenceResultsWriter implements ResultsWriter {
    private static final CSVFormat FORMAT = CSVFormat.TDF.builder().setRecordSeparator("\n")
            .setHeader("CASE_ID", "SAMPLE_PROBAND", "SAMPLE_HPO_IDS", "SCORER", "ENTREZ_GENE_ID", "GENE_SYMBOL", "MODEL_ID", "ORGANISM",
                    "RAW_MODEL_COMPATIBILITY", "RAW_TRANSD_TRIPLE", "NATIVE_PRIORITY_RESULT_SCORE",
                    "MODEL_EVIDENCE_MISSING", "TRANSD_GENE_MISSING", "EMBEDDINGS", "TRANSD_BUNDLE",
                    "EMBEDDINGS_SHA256", "TRANSD_BUNDLE_SHA256")
            .build();

    @Override
    public void writeFile(AnalysisResults results, OutputSettings settings) {
        var out = settings.makeOutputFilePath(results.sample().vcfPath(), OutputFormat.TSV_NATIVE_PHENOTYPE);
        try { Files.writeString(out, writeString(results, settings), StandardCharsets.UTF_8); }
        catch (IOException e) { throw new IllegalStateException("Unable to write native phenotype evidence: " + out, e); }
    }

    @Override
    public String writeString(AnalysisResults results, OutputSettings settings) {
        var options = results.analysis().learnedScoringOptions();
        if (options.method() == org.monarchinitiative.exomiser.core.analysis.LearnedScoringOptions.Method.PHENODIGM)
            throw new IllegalArgumentException("TSV_NATIVE_PHENOTYPE requires learned phenotype scorer");
        TransdCheckpoint checkpoint = null;
        if (options.transdBundle() != null) {
            try { checkpoint = TransdCheckpoint.load(options.transdBundle()); }
            catch (IOException e) { throw new IllegalArgumentException("Unable to load TransD bundle for native provenance", e); }
        }
        StringWriter output = new StringWriter();
        String embeddingsDigest = ScoringResourceDigest.sha256(options.embeddings());
        String bundleDigest = ScoringResourceDigest.transdBundle(options.transdBundle());
        try (CSVPrinter printer = new CSVPrinter(output, FORMAT)) {
            for (var gene : settings.filterGenesForOutput(results.genes())) {
                var hiPhive = gene.getPriorityResult(HiPhivePriorityResult.class);
                var phive = gene.getPriorityResult(PhivePriorityResult.class);
                List<GeneModelPhenotypeMatch> matches = hiPhive != null ? hiPhive.phenotypeEvidence()
                        : phive != null && phive.geneModelPhenotypeMatch() != null ? List.of(phive.geneModelPhenotypeMatch()) : List.of();
                // PHIVE uses a synthetic 0.6 mouse-model match when no database model exists.
                // It affects native ranking but is not raw model compatibility evidence.
                if (hiPhive == null && phive != null && matches.size() == 1 && isPhiveNoMouseModel(matches.getFirst())) {
                    matches = List.of();
                }
                double nativeScore = hiPhive != null ? hiPhive.score() : phive != null ? phive.score() : 0;
                String id = gene.geneIdentifier().entrezId();
                Double triple = checkpoint != null && checkpoint.hasGene(id) && checkpoint.hasCase(options.caseId())
                        ? checkpoint.scoreGeneCase(id, options.caseId()) : null;
                boolean geneMissing = checkpoint != null && !checkpoint.hasGene(id);
                if (matches.isEmpty()) {
                    printer.printRecord(options.caseId(), results.sample().probandSampleName(),
                            String.join(",", results.sample().hpoIds()), options.method(), id, gene.geneSymbol(), "", "", "", triple,
                            nativeScore, true, geneMissing, options.embeddings(), options.transdBundle(), embeddingsDigest, bundleDigest);
                } else {
                    for (var match : matches) {
                        printer.printRecord(options.caseId(), results.sample().probandSampleName(),
                                String.join(",", results.sample().hpoIds()), options.method(), id, gene.geneSymbol(), match.model().id(),
                                match.organism(), match.score(), triple, nativeScore, false,
                                geneMissing, options.embeddings(), options.transdBundle(), embeddingsDigest, bundleDigest);
                    }
                }
            }
        } catch (IOException e) { throw new IllegalStateException("Unable to format native phenotype evidence", e); }
        return output.toString();
    }

    private static boolean isPhiveNoMouseModel(GeneModelPhenotypeMatch match) {
        return match.model() instanceof GeneOrthologModel model
                && model.modelId().isEmpty()
                && model.modelGeneId().isEmpty()
                && model.modelGeneSymbol().isEmpty()
                && model.phenotypeIds().isEmpty()
                && match.bestPhenotypeMatches().isEmpty()
                && match.score() == (double) 0.6f;
    }
}
