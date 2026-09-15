package org.monarchinitiative.exomiser.core.analysis;

import org.junit.jupiter.api.Test;
import org.monarchinitiative.exomiser.api.v1.AnalysisProto;
import static org.junit.jupiter.api.Assertions.*;

class LearnedScoringOptionsTest {
    @Test
    void realJobYamlKeepsLearnedOptionsThroughJobReader() {
        String yaml = """
                sample:
                  genomeAssembly: GRCh38
                  vcf: case.vcf
                  hpoIds: [HP:0001250]
                analysis:
                  phenotypeScorer:
                    method: INDIGENA
                    embeddings: /models/embeddings.tsv
                    caseId: case_1
                  weightedVariants:
                    phenotypeWeight: 0.61
                    genePhenotypes: /data/gene_phenotypes_mp.csv
                outputOptions:
                  outputFormats: [TSV_VARIANT, TSV_TABLE3_WEIGHTED]
                """;
        var job = JobReader.readJob(yaml);
        var options = LearnedScoringOptions.from(job.getAnalysis());
        assertEquals(LearnedScoringOptions.Method.INDIGENA, options.method());
        assertEquals("case_1", options.caseId());
        assertEquals(0.61, options.weighted().phenotypeWeight());
        assertTrue(job.getOutputOptions().getOutputFormatsList().contains("TSV_TABLE3_WEIGHTED"));
    }

    @Test
    void defaultsPreservePhenodigmAndInvalidCombinationsFail() {
        assertEquals(LearnedScoringOptions.Method.PHENODIGM,
                LearnedScoringOptions.from(AnalysisProto.Analysis.getDefaultInstance()).method());
        var unknown = AnalysisProto.Analysis.newBuilder().setPhenotypeScorer(
                AnalysisProto.PhenotypeScorerOptions.newBuilder().setMethod("EMBEDPVP")).build();
        assertThrows(IllegalArgumentException.class, () -> LearnedScoringOptions.from(unknown));
        var incomplete = AnalysisProto.Analysis.newBuilder().setPhenotypeScorer(
                AnalysisProto.PhenotypeScorerOptions.newBuilder().setMethod("EMBEDPVP_TRANSD")
                        .setTransdBundle("/bundle")).build();
        assertThrows(IllegalArgumentException.class, () -> LearnedScoringOptions.from(incomplete));
    }
}
