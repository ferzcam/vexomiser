package org.monarchinitiative.exomiser.core.writers;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;
import org.monarchinitiative.exomiser.core.analysis.Analysis;
import org.monarchinitiative.exomiser.core.analysis.AnalysisResults;
import org.monarchinitiative.exomiser.core.analysis.LearnedScoringOptions;
import org.monarchinitiative.exomiser.core.model.Gene;

import java.nio.file.Files;
import java.nio.file.Path;
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

    private static String field(String headerLine, String rowLine, String name) {
        String[] header = headerLine.split("\\t", -1);
        String[] row = rowLine.split("\\t", -1);
        for (int i = 0; i < header.length; i++) if (header[i].equals(name)) return row[i];
        throw new IllegalArgumentException("Missing column: " + name);
    }
}
