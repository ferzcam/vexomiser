package org.monarchinitiative.exomiser.core.phenotype;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.Set;
import static org.junit.jupiter.api.Assertions.*;

class PooledGenePhenotypeScorerTest {
    @TempDir Path temporary;

    @Test
    void scoresPooledAnnotationsAndKeepsExplicitZeroForUnannotatedGene() throws Exception {
        Path emb = temporary.resolve("emb.tsv");
        Files.writeString(emb, "http://purl.obolibrary.org/obo/HP_1\t1,0\n" +
                "http://purl.obolibrary.org/obo/MP_1\t1,0\n" +
                "http://purl.obolibrary.org/obo/MP_2\t0,1\n");
        Path csv = temporary.resolve("gene.csv");
        Files.writeString(csv, "gene,phenotype\nhttp://mowl.borg/1,http://purl.obolibrary.org/obo/MP_1\n" +
                "http://mowl.borg/1,http://purl.obolibrary.org/obo/MP_2\n");
        var scorer = new PooledGenePhenotypeScorer(IndigenaEmbeddings.load(emb), csv);
        var scores = scorer.score(Set.of("1", "2"), List.of("HP:1"));
        assertEquals(0.673294, scores.get("1"), 1e-6);
        assertEquals(0.0, scores.get("2"));
    }

    @Test
    void sixDecimalRoundingUsesExactBinaryValue() {
        assertEquals(0.000003, PooledGenePhenotypeScorer.sixDecimal(0.0000025));
    }
}
