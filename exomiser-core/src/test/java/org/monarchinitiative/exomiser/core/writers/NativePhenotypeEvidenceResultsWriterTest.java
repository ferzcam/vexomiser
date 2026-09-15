package org.monarchinitiative.exomiser.core.writers;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;
import org.monarchinitiative.exomiser.core.analysis.Analysis;
import org.monarchinitiative.exomiser.core.analysis.AnalysisResults;
import org.monarchinitiative.exomiser.core.analysis.LearnedScoringOptions;
import org.monarchinitiative.exomiser.core.model.Gene;
import org.monarchinitiative.exomiser.core.phenotype.Organism;
import org.monarchinitiative.exomiser.core.prioritisers.PhivePriorityResult;
import org.monarchinitiative.exomiser.core.prioritisers.model.GeneModelPhenotypeMatch;
import org.monarchinitiative.exomiser.core.prioritisers.model.GeneOrthologModel;

import java.nio.file.Files;
import java.nio.file.Path;
import java.security.MessageDigest;
import java.util.List;
import static org.junit.jupiter.api.Assertions.*;

class NativePhenotypeEvidenceResultsWriterTest {
    @TempDir Path temporary;

    @Test
    void geneWithoutDatabaseModelProducesExplicitMissingEvidenceRow() throws Exception {
        Path embeddings = temporary.resolve("emb.tsv");
        Files.writeString(embeddings, "http://purl.obolibrary.org/obo/HP_1\t1,0\n");
        var options = new LearnedScoringOptions(LearnedScoringOptions.Method.INDIGENA,
                embeddings, null, "case_1", null);
        var results = AnalysisResults.builder()
                .analysis(Analysis.builder().learnedScoringOptions(options).build())
                .genes(List.of(new Gene("UNMODELED", 123456))).build();
        String output = new NativePhenotypeEvidenceResultsWriter().writeString(results, OutputSettings.defaults());
        assertEquals(2, output.lines().count());
        String header = output.lines().findFirst().orElseThrow();
        String row = output.lines().skip(1).findFirst().orElseThrow();
        assertEquals("true", field(header, row, "MODEL_EVIDENCE_MISSING"));
        assertEquals("", field(header, row, "RAW_MODEL_COMPATIBILITY"));
        assertEquals("case_1", field(header, row, "CASE_ID"));
    }

    @Test
    void phiveFallbackKeepsNativeScoreAndTripleButHasNoRawModelEvidence() throws Exception {
        Path bundle = transdBundle();
        var options = new LearnedScoringOptions(LearnedScoringOptions.Method.EMBEDPVP_TRANSD,
                null, bundle, "case_1", null);
        Gene gene = new Gene("UNMODELED", 123456);
        GeneOrthologModel placeholder = new GeneOrthologModel("", Organism.MOUSE, 123456,
                "UNMODELED", "", "", List.of());
        double fallback = (double) 0.6f;
        gene.addPriorityResult(new PhivePriorityResult(123456, "UNMODELED", fallback,
                new GeneModelPhenotypeMatch(fallback, placeholder, List.of())));
        String output = new NativePhenotypeEvidenceResultsWriter().writeString(
                results(options, gene), OutputSettings.defaults());
        String header = output.lines().findFirst().orElseThrow();
        String row = output.lines().skip(1).findFirst().orElseThrow();
        assertEquals(2, output.lines().count());
        assertEquals("true", field(header, row, "MODEL_EVIDENCE_MISSING"));
        assertEquals("", field(header, row, "RAW_MODEL_COMPATIBILITY"));
        assertEquals("", field(header, row, "MODEL_ID"));
        assertEquals(Double.toString(fallback), field(header, row, "NATIVE_PRIORITY_RESULT_SCORE"));
        assertEquals("-1.0", field(header, row, "RAW_TRANSD_TRIPLE"));
        assertEquals("false", field(header, row, "TRANSD_GENE_MISSING"));
    }

    @Test
    void phiveRealMouseModelExportsRawCompatibility() throws Exception {
        Path bundle = transdBundle();
        var options = new LearnedScoringOptions(LearnedScoringOptions.Method.EMBEDPVP_TRANSD,
                null, bundle, "case_1", null);
        Gene gene = new Gene("MODELED", 123457);
        GeneOrthologModel model = new GeneOrthologModel("MGI:123", Organism.MOUSE, 123457,
                "MODELED", "MGI:123", "MouseGene", List.of("MP:0001"));
        gene.addPriorityResult(new PhivePriorityResult(123457, "MODELED", 0.72,
                new GeneModelPhenotypeMatch(0.72, model, List.of())));
        String output = new NativePhenotypeEvidenceResultsWriter().writeString(
                results(options, gene), OutputSettings.defaults());
        String header = output.lines().findFirst().orElseThrow();
        String row = output.lines().skip(1).findFirst().orElseThrow();
        assertEquals("false", field(header, row, "MODEL_EVIDENCE_MISSING"));
        assertEquals("0.72", field(header, row, "RAW_MODEL_COMPATIBILITY"));
        assertEquals("MGI:123", field(header, row, "MODEL_ID"));
        assertEquals("0.72", field(header, row, "NATIVE_PRIORITY_RESULT_SCORE"));
        assertEquals("-0.25", field(header, row, "RAW_TRANSD_TRIPLE"));
    }

    private static AnalysisResults results(LearnedScoringOptions options, Gene gene) {
        return AnalysisResults.builder().analysis(Analysis.builder().learnedScoringOptions(options).build())
                .genes(List.of(gene)).build();
    }

    private Path transdBundle() throws Exception {
        Path bundle = temporary.resolve("transd");
        Files.createDirectories(bundle);
        Files.writeString(bundle.resolve("entities.tsv"),
                "http://mowl.borg/123456\t1.0,0.0\t0.0,0.0\n" +
                "http://mowl.borg/123457\t0.5,0.0\t0.0,0.0\n" +
                "http://mowl.borg/case_1\t0.0,0.0\t0.0,0.0\n");
        Files.writeString(bundle.resolve("relations.tsv"),
                "http://mowl.borg/associated_with\t0.0,0.0\t0.0,0.0\n");
        Files.writeString(bundle.resolve("manifest.properties"),
                "format=pykeen-transd-v1\npykeen_version=1.11.1\n" +
                "model_state_sha256=" + "0".repeat(64) + "\n" +
                "entity_dimension=2\nrelation_dimension=2\nentity_count=3\nrelation_count=1\n" +
                "entities_sha256=" + sha256(bundle.resolve("entities.tsv")) + "\n" +
                "relations_sha256=" + sha256(bundle.resolve("relations.tsv")) + "\n" +
                "norm_p=2\npower_norm=true\nprojection_maxnorm=1\ninverse_triples=false\n" +
                "predict_with_sigmoid=false\nassociation_relation=http://mowl.borg/associated_with\n");
        return bundle;
    }

    private static String sha256(Path file) throws Exception {
        byte[] digest = MessageDigest.getInstance("SHA-256").digest(Files.readAllBytes(file));
        StringBuilder result = new StringBuilder(64);
        for (byte value : digest) result.append(String.format("%02x", value & 0xff));
        return result.toString();
    }

    private static String field(String headerLine, String rowLine, String name) {
        String[] header = headerLine.split("\\t", -1);
        String[] row = rowLine.split("\\t", -1);
        for (int i = 0; i < header.length; i++) if (header[i].equals(name)) return row[i];
        throw new IllegalArgumentException("Missing column: " + name);
    }
}
