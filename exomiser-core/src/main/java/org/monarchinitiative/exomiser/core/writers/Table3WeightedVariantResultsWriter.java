package org.monarchinitiative.exomiser.core.writers;

import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVPrinter;
import org.monarchinitiative.exomiser.core.analysis.AnalysisResults;
import org.monarchinitiative.exomiser.core.analysis.LearnedScoringOptions;
import org.monarchinitiative.exomiser.core.analysis.score.WeightedVariantScorer;
import org.monarchinitiative.exomiser.core.phenotype.IndigenaEmbeddings;
import org.monarchinitiative.exomiser.core.phenotype.PooledGenePhenotypeScorer;
import org.monarchinitiative.exomiser.core.phenotype.TransdCheckpoint;
import org.monarchinitiative.exomiser.core.phenotype.ScoringResourceDigest;
import org.monarchinitiative.exomiser.core.model.pathogenicity.PathogenicityScore;

import java.io.IOException;
import java.io.StringWriter;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.math.BigDecimal;
import java.math.RoundingMode;

/** Table 3's distinct variant ranks over exactly the exported native TSV candidate membership. */
public final class Table3WeightedVariantResultsWriter implements ResultsWriter {
    private static final CSVFormat FORMAT = CSVFormat.TDF.builder()
            .setRecordSeparator("\n")
            .setHeader("CASE_ID", "SAMPLE_PROBAND", "SAMPLE_HPO_IDS", "SCORER", "CANDIDATE_SOURCE", "CONTIG", "START", "REF", "ALT",
                    "ENTREZ_GENE_IDS", "RAW_MS", "RAW_GP_MAX_PATH", "MS_MISSING", "GP_MISSING",
                    "IMPUTED_MS", "IMPUTED_GP", "NORMALIZED_MS", "NORMALIZED_GP", "WEIGHT",
                    "COMBINED_SCORE", "VARIANT_RANK_MID", "VARIANT_RANK_STRICT",
                    "EMBEDDINGS", "TRANSD_BUNDLE", "GENE_PHENOTYPES", "PATIENT_HPO",
                    "MS_MISSING_POLICY", "GP_MISSING_POLICY", "EMBEDDINGS_SHA256", "TRANSD_BUNDLE_SHA256",
                    "GENE_PHENOTYPES_SHA256", "PATIENT_HPO_SHA256", "PATHOGENICITY_SOURCES")
            .build();

    @Override
    public void writeFile(AnalysisResults results, OutputSettings settings) {
        Path out = settings.makeOutputFilePath(results.sample().vcfPath(), OutputFormat.TSV_TABLE3_WEIGHTED);
        try {
            Files.writeString(out, writeString(results, settings), StandardCharsets.UTF_8);
        } catch (IOException e) {
            throw new IllegalStateException("Unable to write weighted variant output: " + out, e);
        }
    }

    @Override
    public String writeString(AnalysisResults results, OutputSettings settings) {
        LearnedScoringOptions options = results.analysis().learnedScoringOptions();
        if (options.weighted() == null)
            throw new IllegalArgumentException("TSV_TABLE3_WEIGHTED requires analysis.weightedVariants");
        List<WeightedVariantScorer.VariantInput> rows = candidateRows(results, settings);
        Set<String> genes = new HashSet<>();
        rows.forEach(row -> genes.addAll(row.geneIds()));
        Map<String, Double> scores = phenotypeScores(options, genes, results);
        var weighted = options.weighted();
        var scoreSettings = new WeightedVariantScorer.Settings(weighted.phenotypeWeight(),
                weighted.pathogenicityMissing(), weighted.phenotypeMissing());
        List<WeightedVariantScorer.VariantScore> ranks = new WeightedVariantScorer().score(rows, scores, scoreSettings);
        String embeddingsDigest = ScoringResourceDigest.sha256(options.embeddings());
        String transdDigest = ScoringResourceDigest.transdBundle(options.transdBundle());
        String annotationDigest = ScoringResourceDigest.sha256(weighted.genePhenotypes());
        String patientDigest = ScoringResourceDigest.sha256(weighted.patientHpo());
        StringWriter output = new StringWriter();
        try (CSVPrinter printer = new CSVPrinter(output, FORMAT)) {
            for (var rank : ranks) {
                var key = rank.key();
                printer.printRecord(options.caseId(), results.sample().probandSampleName(),
                        String.join(",", results.sample().hpoIds()), options.method(),
                        "GeneScoreRanker.rankedVariants TSV_VARIANT",
                        key.contig(), key.start(), key.ref(), key.alt(),
                        String.join(",", rank.geneIds().stream().sorted().toList()),
                        rank.rawPhenotypeScore(), rank.rawPathogenicityScore(),
                        rank.phenotypeMissing(), rank.pathogenicityMissing(),
                        rank.imputedPhenotypeScore(), rank.imputedPathogenicityScore(),
                        rank.normalizedPhenotypeScore(), rank.normalizedPathogenicityScore(),
                        weighted.phenotypeWeight(), rank.combinedScore(), rank.rank(), rank.rankStrict(),
                        options.embeddings(), options.transdBundle(), weighted.genePhenotypes(), weighted.patientHpo(),
                        weighted.phenotypeMissing(), weighted.pathogenicityMissing(),
                        embeddingsDigest, transdDigest, annotationDigest, patientDigest,
                        String.join(",", results.analysis().pathogenicitySources().stream().map(Enum::name).sorted().toList()));
            }
        } catch (IOException e) {
            throw new IllegalStateException("Unable to format weighted variant output", e);
        }
        return output.toString();
    }

    private static List<WeightedVariantScorer.VariantInput> candidateRows(AnalysisResults results, OutputSettings settings) {
        List<WeightedVariantScorer.VariantInput> rows = new ArrayList<>();
        new GeneScoreRanker(results, settings).rankedVariants().forEach(row -> {
            rows.add(legacyTsvRowInput(row));
        });
        return rows;
    }

    /** Match the Float.toString handoff of the legacy TSV before Python parses MAX_PATH as double. */
    static WeightedVariantScorer.VariantInput legacyTsvRowInput(GeneScoreRanker.RankedVariant row) {
        var ve = row.variantEvaluation();
        var id = row.geneScore().geneIdentifier().entrezId();
        PathogenicityScore maxPath = ve.pathogenicityData().mostPathogenicScore();
        // A null predictor is exported as an empty field; the combiner treats empty or '.' as missing.
        Double gp = maxPath == null ? null : Double.valueOf(Float.toString(maxPath.score()));
        return new WeightedVariantScorer.VariantInput(
                new WeightedVariantScorer.VariantKey(ve.contigName(), ve.start(), ve.ref(), ve.alt()),
                gp, id.isBlank() ? Set.of() : Set.of(id));
    }

    private static Map<String, Double> phenotypeScores(LearnedScoringOptions options, Set<String> genes,
                                                         AnalysisResults results) {
        try {
            return switch (options.method()) {
                case INDIGENA -> {
                    var embeddings = IndigenaEmbeddings.load(options.embeddings());
                    var pooled = new PooledGenePhenotypeScorer(embeddings, options.weighted().genePhenotypes());
                    List<String> terms = patientTerms(options.weighted().patientHpo(), results);
                    yield pooled.score(genes, terms);
                }
                case EMBEDPVP_TRANSD -> {
                    var checkpoint = TransdCheckpoint.load(options.transdBundle());
                    if (!checkpoint.hasCase(options.caseId()))
                        throw new IllegalArgumentException("TransD bundle has no case entity: " + options.caseId());
                    Map<String, Double> raw = new HashMap<>();
                    for (String gene : genes) {
                        if (checkpoint.hasGene(gene)) raw.put(gene, sixDecimal(checkpoint.scoreGeneCase(gene, options.caseId())));
                    }
                    yield raw;
                }
                case PHENODIGM -> throw new IllegalArgumentException("Weighted output requires learned scorer");
            };
        } catch (IOException e) {
            throw new IllegalArgumentException("Unable to load weighted phenotype evidence", e);
        }
    }

    private static List<String> patientTerms(Path controlled, AnalysisResults results) throws IOException {
        if (controlled == null) return results.sample().hpoIds();
        List<String> terms = Files.readAllLines(controlled, StandardCharsets.UTF_8).stream()
                .map(String::trim).filter(line -> !line.isBlank() && !line.startsWith("#")).toList();
        if (terms.isEmpty()) throw new IllegalArgumentException("Controlled patient HPO list is empty: " + controlled);
        return terms;
    }

    /** The historical per-gene phenotype TSV emits f"{score:.6f}" before the combiner parses it. */
    static double sixDecimal(double score) {
        if (!Double.isFinite(score)) throw new IllegalArgumentException("Non-finite phenotype score");
        return new BigDecimal(score).setScale(6, RoundingMode.HALF_EVEN).doubleValue();
    }
}
