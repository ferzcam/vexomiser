package org.monarchinitiative.exomiser.core.analysis;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;
import org.monarchinitiative.exomiser.core.genome.GenomeAnalysisServiceProvider;
import org.monarchinitiative.exomiser.core.genome.TestFactory;
import org.monarchinitiative.exomiser.core.model.Gene;
import org.monarchinitiative.exomiser.core.phenotype.IndigenaModelScorerFactory;
import org.monarchinitiative.exomiser.core.prioritisers.HiPhivePriority;
import org.monarchinitiative.exomiser.core.prioritisers.PriorityFactoryImpl;
import org.monarchinitiative.exomiser.core.prioritisers.service.TestPriorityServiceFactory;
import org.monarchinitiative.exomiser.core.prioritisers.util.DataMatrix;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import static org.junit.jupiter.api.Assertions.*;

class JobParserLearnedInvocationTest {
    @TempDir Path temporary;

    @Test
    void yamlSelectorReachesNativePrioritiserInvocation() throws Exception {
        Path embeddings = temporary.resolve("embeddings.tsv");
        Files.writeString(embeddings, "http://purl.obolibrary.org/obo/HP_0001156\t1,0\n" +
                "http://purl.obolibrary.org/obo/MP_0000031\t1,0\n");
        String yaml = """
                sample:
                  genomeAssembly: GRCh37
                  vcf: case.vcf
                  hpoIds: [HP:0001156]
                analysis:
                  phenotypeScorer:
                    method: INDIGENA
                    embeddings: %s
                  steps:
                    - hiPhivePrioritiser: {runParams: mouse}
                """.formatted(embeddings);
        var job = JobReader.readJob(yaml);
        var service = TestPriorityServiceFactory.testPriorityService();
        var parser = new JobParser(new GenomeAnalysisServiceProvider(TestFactory.buildDefaultHg19GenomeAnalysisService()),
                new PriorityFactoryImpl(service, DataMatrix.empty(), temporary),
                TestPriorityServiceFactory.testOntologyService());
        Analysis analysis = parser.parseAnalysis(job);
        assertEquals(LearnedScoringOptions.Method.INDIGENA, analysis.learnedScoringOptions().method());
        HiPhivePriority prioritiser = assertInstanceOf(HiPhivePriority.class, analysis.mainPrioritiser());
        var factoryField = HiPhivePriority.class.getDeclaredField("modelScorerFactory");
        factoryField.setAccessible(true);
        assertInstanceOf(IndigenaModelScorerFactory.class, factoryField.get(prioritiser));
        var result = prioritiser.prioritise(List.of("HP:0001156"), List.of(new Gene("FGFR2", 2263)))
                .findFirst().orElseThrow();
        // MP:0000031 is present in the test FGFR2 mouse model. The sigmoid(dot=1) learned
        // match survives the CLI JobProto -> JobParser -> native HiPhive invocation path.
        assertEquals(1.0 / (1.0 + Math.exp(-1)), result.mouseScore(), 1e-6);
    }

    @Test
    void learnedSelectorWithoutCompatiblePrioritiserFailsInsteadOfBeingIgnored() {
        String yaml = """
                sample:
                  genomeAssembly: GRCh37
                  vcf: case.vcf
                  hpoIds: [HP:0001156]
                analysis:
                  phenotypeScorer:
                    method: INDIGENA
                    embeddings: /models/emb.tsv
                  steps:
                    - phenixPrioritiser: {}
                """;
        var service = TestPriorityServiceFactory.testPriorityService();
        var parser = new JobParser(new GenomeAnalysisServiceProvider(TestFactory.buildDefaultHg19GenomeAnalysisService()),
                new PriorityFactoryImpl(service, DataMatrix.empty(), temporary),
                TestPriorityServiceFactory.testOntologyService());
        var exception = assertThrows(IllegalArgumentException.class,
                () -> parser.parseAnalysis(JobReader.readJob(yaml)));
        assertTrue(exception.getMessage().contains("hiPhivePrioritiser or phivePrioritiser"));
    }
}
