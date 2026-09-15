package org.monarchinitiative.exomiser.core.analysis.score;

import org.junit.jupiter.api.Test;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import static org.junit.jupiter.api.Assertions.*;
import static org.monarchinitiative.exomiser.core.analysis.score.WeightedVariantScorer.*;

class WeightedVariantScorerTest {

    private final WeightedVariantScorer scorer = new WeightedVariantScorer();

    private static VariantInput row(String contig, int position, Double gp, String... genes) {
        return new VariantInput(new VariantKey(contig, position, "A", "G"), gp, Set.of(genes));
    }

    private static VariantScore at(List<VariantScore> scores, int position) {
        return scores.stream().filter(s -> s.key().start() == position).findFirst().orElseThrow();
    }

    @Test
    void duplicateRowsUnionGenesAndKeepMaximumAvailablePathogenicity() {
        List<VariantScore> scores = scorer.score(List.of(
                row("chr1", 1, null, "negative"),
                row("1", 1, -0.2, "positive"),
                row("1", 1, 0.8, "other"),
                row("1", 2, 0.0, "negative")),
                Map.of("negative", -5.0, "positive", -2.0, "other", -3.0), Settings.table3(0.6));

        assertEquals(2, scores.size());
        VariantScore first = at(scores, 1);
        assertEquals(Set.of("negative", "positive", "other"), first.geneIds());
        assertEquals(-2.0, first.rawPhenotypeScore());
        assertEquals(0.8, first.rawPathogenicityScore());
        assertFalse(first.phenotypeMissing());
        assertFalse(first.pathogenicityMissing());
        assertEquals(1.0, first.normalizedPhenotypeScore());
        assertEquals(1.0, first.normalizedPathogenicityScore());
        assertEquals(1.0, first.rank());
    }

    @Test
    void explicitZeroAndAllMissingRowsStayDistinct() {
        List<VariantScore> scores = scorer.score(List.of(
                row("1", 1, 0.0, "zero"), row("1", 2, null, "missing"),
                row("1", 3, 1.0, "high")),
                Map.of("zero", 0.0, "high", 4.0), Settings.table3(0.5));

        assertEquals(3, scores.size());
        assertEquals(0.0, at(scores, 1).rawPhenotypeScore());
        assertFalse(at(scores, 1).phenotypeMissing());
        assertFalse(at(scores, 1).pathogenicityMissing());
        assertNull(at(scores, 2).rawPhenotypeScore());
        assertNull(at(scores, 2).rawPathogenicityScore());
        assertTrue(at(scores, 2).phenotypeMissing());
        assertTrue(at(scores, 2).pathogenicityMissing());
        assertEquals(2.0, at(scores, 2).imputedPhenotypeScore());
        assertEquals(0.0, at(scores, 2).imputedPathogenicityScore());
        assertEquals(0.5, at(scores, 2).normalizedPhenotypeScore());
    }

    @Test
    void mediansCountScoreableVariantsRatherThanDistinctGenesOrInputRows() {
        List<VariantScore> scores = scorer.score(List.of(
                row("1", 1, 1.0, "shared"), row("1", 1, 1.0, "shared"),
                row("1", 2, 3.0, "shared"), row("1", 3, 8.0, "other"),
                row("1", 4, null, "absent")),
                Map.of("shared", -4.0, "other", -1.0),
                new Settings(0.5, MissingPolicy.MEDIAN, MissingPolicy.MEDIAN));

        assertEquals(4, scores.size());
        assertEquals(3.0, at(scores, 4).imputedPathogenicityScore());
        assertEquals(-4.0, at(scores, 4).imputedPhenotypeScore());
        assertEquals(-4.0, at(scores, 2).rawPhenotypeScore());
    }

    @Test
    void constantColumnsAndExactTiesHaveAverageAndStrictRanks() {
        List<VariantScore> scores = scorer.score(List.of(
                row("1", 3, 2.0, "g"), row("1", 1, 2.0, "g"),
                row("1", 2, 2.0, "g")), Map.of("g", -3.0), Settings.table3(0.6));

        assertEquals(List.of(1, 2, 3), scores.stream().map(s -> s.key().start()).toList());
        for (VariantScore score : scores) {
            assertEquals(0.0, score.normalizedPhenotypeScore());
            assertEquals(0.0, score.normalizedPathogenicityScore());
            assertEquals(0.0, score.combinedScore());
            assertEquals(2.0, score.rank());
            assertEquals(1.0, score.rankStrict());
        }
    }

    @Test
    void endpointWeightsRankOnlyTheirRespectiveComponent() {
        List<VariantInput> rows = List.of(row("1", 1, 0.0, "best"), row("1", 2, 1.0, "worst"));
        Map<String, Double> genes = Map.of("best", 10.0, "worst", -10.0);
        assertEquals(1, scorer.score(rows, genes, Settings.table3(1.0)).get(0).key().start());
        assertEquals(2, scorer.score(rows, genes, Settings.table3(0.0)).get(0).key().start());
    }

    @Test
    void rejectsInvalidNumbersAndSettings() {
        assertThrows(IllegalArgumentException.class, () -> row("1", 1, Double.NaN, "g"));
        assertThrows(IllegalArgumentException.class, () -> row("1", 1, Double.POSITIVE_INFINITY, "g"));
        assertThrows(IllegalArgumentException.class, () -> Settings.table3(-0.1));
        assertThrows(IllegalArgumentException.class, () -> Settings.table3(Double.NaN));
        assertThrows(NullPointerException.class, () -> new Settings(0.5, null, MissingPolicy.ZERO));
        assertThrows(IllegalArgumentException.class, () -> scorer.score(List.of(row("1", 1, 0.0, "g")),
                Map.of("g", Double.NEGATIVE_INFINITY), Settings.table3(0.5)));
        assertThrows(IllegalArgumentException.class, () -> scorer.score(List.of(
                        row("1", 1, -Double.MAX_VALUE, "g"), row("1", 2, Double.MAX_VALUE, "g")),
                Map.of("g", 0.0), Settings.table3(0.5)));
    }

    @Test
    void inputAndOutputCollectionsAreDefensivelyImmutable() {
        Set<String> mutableGenes = new HashSet<>(Set.of("g"));
        VariantInput input = new VariantInput(new VariantKey("chr1", 1, "A", "G"), null, mutableGenes);
        mutableGenes.add("later");
        Map<String, Double> table = new HashMap<>();
        table.put("g", 1.0);
        List<VariantScore> scores = scorer.score(new ArrayList<>(List.of(input)), table, Settings.table3(0.6));
        table.put("g", 9.0);
        assertEquals(Set.of("g"), input.geneIds());
        assertEquals(Set.of("g"), scores.get(0).geneIds());
        assertEquals(1.0, scores.get(0).rawPhenotypeScore());
        assertThrows(UnsupportedOperationException.class, () -> input.geneIds().add("x"));
        assertThrows(UnsupportedOperationException.class, () -> scores.get(0).geneIds().add("x"));
        assertThrows(UnsupportedOperationException.class, () -> scores.add(scores.get(0)));
    }
}
