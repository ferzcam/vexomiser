package org.monarchinitiative.exomiser.core.writers;

import de.charite.compbio.jannovar.mendel.ModeOfInheritance;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;
import org.monarchinitiative.exomiser.core.analysis.Analysis;
import org.monarchinitiative.exomiser.core.analysis.AnalysisResults;
import org.monarchinitiative.exomiser.core.analysis.LearnedScoringOptions;
import org.monarchinitiative.exomiser.core.analysis.score.WeightedVariantScorer;
import org.monarchinitiative.exomiser.core.filters.FilterResult;
import org.monarchinitiative.exomiser.core.filters.FilterType;
import org.monarchinitiative.exomiser.core.model.Gene;
import org.monarchinitiative.exomiser.core.model.GeneScore;
import org.monarchinitiative.exomiser.core.model.SampleGenotype;
import org.monarchinitiative.exomiser.core.genome.TestFactory;
import org.monarchinitiative.exomiser.core.genome.TestVariantFactory;
import org.monarchinitiative.exomiser.core.model.VariantEvaluation;
import org.monarchinitiative.exomiser.core.model.pathogenicity.PathogenicityData;
import org.monarchinitiative.exomiser.core.model.pathogenicity.PathogenicityScore;
import org.monarchinitiative.exomiser.core.model.pathogenicity.PathogenicitySource;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.EnumSet;
import java.util.List;
import static org.junit.jupiter.api.Assertions.*;

class Table3WeightedVariantResultsWriterTest {
    @TempDir Path temporary;

    @Test
    void candidateRowsMatchLegacyTsvAndVariantsInOneGeneGetIndependentRanks() throws Exception {
        Path emb = temporary.resolve("emb.tsv");
        Files.writeString(emb, "http://purl.obolibrary.org/obo/HP_1\t1,0\n" +
                "http://purl.obolibrary.org/obo/MP_1\t1,0\n");
        Path csv = temporary.resolve("gene.csv");
        Files.writeString(csv, "gene,phenotype\nhttp://mowl.borg/2263,http://purl.obolibrary.org/obo/MP_1\n");
        LearnedScoringOptions options = new LearnedScoringOptions(LearnedScoringOptions.Method.INDIGENA,
                emb, null, "case_1", new LearnedScoringOptions.Weighted(0.61, csv, null,
                WeightedVariantScorer.MissingPolicy.ZERO, WeightedVariantScorer.MissingPolicy.MEDIAN));
        Gene gene = TestFactory.newGeneFGFR2();
        VariantEvaluation first = variant(123256214, 0.2f);
        VariantEvaluation second = variant(123256215, 0.9f);
        VariantEvaluation missing = variant(123256216, null);
        gene.addVariant(first);
        gene.addVariant(second);
        gene.addVariant(missing);
        gene.addGeneScore(GeneScore.builder().geneIdentifier(gene.geneIdentifier())
                .combinedScore(1.0).phenotypeScore(0.7).variantScore(0.5)
                .contributingVariants(List.of(first, second, missing))
                .modeOfInheritance(ModeOfInheritance.AUTOSOMAL_DOMINANT).build());
        gene.setCompatibleInheritanceModes(EnumSet.of(ModeOfInheritance.AUTOSOMAL_DOMINANT));
        Gene otherGene = TestFactory.newGeneSHH();
        VariantEvaluation duplicate = variant(123256214, 0.3f);
        otherGene.addVariant(duplicate);
        otherGene.addGeneScore(GeneScore.builder().geneIdentifier(otherGene.geneIdentifier())
                .combinedScore(1.0).phenotypeScore(0.7).variantScore(0.5)
                .contributingVariants(List.of(duplicate))
                .modeOfInheritance(ModeOfInheritance.AUTOSOMAL_DOMINANT).build());
        otherGene.setCompatibleInheritanceModes(EnumSet.of(ModeOfInheritance.AUTOSOMAL_DOMINANT));
        Analysis analysis = Analysis.builder().learnedScoringOptions(options).build();
        AnalysisResults results = AnalysisResults.builder().analysis(analysis).genes(List.of(gene, otherGene)).build();
        OutputSettings settings = OutputSettings.defaults();
        String legacy = new TsvVariantResultsWriter().writeString(results, settings);
        String weighted = new Table3WeightedVariantResultsWriter().writeString(results, settings);
        List<String> golden = Files.readAllLines(Path.of("src/test/resources/weighted/golden-candidate-rows.tsv"))
                .stream().skip(1).sorted().toList();
        List<String> actualLegacy = legacy.lines().skip(1).map(row -> String.join("\t",
                field(legacy, row, "CONTIG"), field(legacy, row, "START"),
                field(legacy, row, "REF"), field(legacy, row, "ALT"),
                field(legacy, row, "ENTREZ_GENE_ID"), field(legacy, row, "MAX_PATH"))).sorted().toList();
        assertEquals(golden, actualLegacy);
        assertEquals(5, legacy.lines().count()); // duplicate assignment retains a legacy TSV row
        assertEquals(4, weighted.lines().count()); // collapsed to three genomic variants
        String[] ranked = weighted.lines().skip(1).toArray(String[]::new);
        assertEquals(3, ranked.length);
        assertTrue(ranked[0].contains("123256215"));
        assertTrue(ranked[1].contains("123256214"));
        assertTrue(ranked[2].contains("123256216"));
        assertTrue(ranked[1].contains("2263") && ranked[1].contains(otherGene.geneIdentifier().entrezId()));
        assertEquals(0.3, Double.parseDouble(field(weighted, ranked[1], "RAW_GP_MAX_PATH")));
        assertEquals("true", field(weighted, ranked[2], "GP_MISSING"));
        assertEquals("", field(weighted, ranked[2], "RAW_GP_MAX_PATH"));
        assertTrue(ranked[0].contains("\t1.0\t"));
        assertTrue(ranked[1].contains("\t2.0\t"));
        assertTrue(ranked[2].contains("\t3.0\t"));
        String[] legacyRows = legacy.lines().skip(1).toArray(String[]::new);
        String legacyGp = java.util.Arrays.stream(legacyRows).filter(line -> line.contains("123256214") &&
                line.contains(otherGene.geneIdentifier().entrezId())).findFirst()
                .map(line -> field(legacy, line, "MAX_PATH")).orElseThrow();
        assertEquals(Double.parseDouble(legacyGp), Double.parseDouble(field(weighted, ranked[1], "RAW_GP_MAX_PATH")));
    }

    @Test
    void nativePhenotypeDependentGeneScoreCanChangeExportedCandidateMembership() {
        OutputSettings settings = OutputSettings.builder().outputContributingVariantsOnly(true).build();
        AnalysisResults low = resultsWithNativeScore(0.0);
        AnalysisResults high = resultsWithNativeScore(1.0);
        assertEquals(0, new GeneScoreRanker(low, settings).rankedVariants().count());
        assertEquals(1, new GeneScoreRanker(high, settings).rankedVariants().count());
    }

    @Test
    void phenotypeTablePrecisionAppliesToNegativeTransdScores() {
        assertEquals(-0.123457, Table3WeightedVariantResultsWriter.sixDecimal(-0.1234567));
        assertEquals(0.123457, Table3WeightedVariantResultsWriter.sixDecimal(0.1234567));
        // 0.0000025 is binary-above the exact halfway point. Python's .6f rounds the
        // binary value to 0.000003; BigDecimal.valueOf would round the shortest decimal
        // representation's even tie to 0.000002.
        assertEquals(0.000003, Table3WeightedVariantResultsWriter.sixDecimal(0.0000025));
    }

    private static AnalysisResults resultsWithNativeScore(double score) {
        Gene gene = TestFactory.newGeneFGFR2();
        VariantEvaluation ve = variant(123256220, 0.5f);
        gene.addVariant(ve);
        gene.addGeneScore(GeneScore.builder().geneIdentifier(gene.geneIdentifier())
                .combinedScore(score).phenotypeScore(score).variantScore(0.5)
                .contributingVariants(List.of(ve))
                .modeOfInheritance(ModeOfInheritance.AUTOSOMAL_DOMINANT).build());
        gene.setCompatibleInheritanceModes(EnumSet.of(ModeOfInheritance.AUTOSOMAL_DOMINANT));
        return AnalysisResults.builder().genes(List.of(gene)).build();
    }

    private static String field(String tsv, String row, String name) {
        String[] header = tsv.lines().findFirst().orElseThrow().split("\\t", -1);
        String[] cells = row.split("\\t", -1);
        for (int i = 0; i < header.length; i++) if (header[i].equals(name)) return cells[i];
        throw new IllegalArgumentException("Missing column: " + name);
    }

    private static VariantEvaluation variant(int start, Float predictor) {
        VariantEvaluation ve = TestVariantFactory.buildVariant(10, start, "A", "G", SampleGenotype.het(), 30, 2.2);
        ve.addFilterResult(FilterResult.pass(FilterType.VARIANT_EFFECT_FILTER));
        ve.setCompatibleInheritanceModes(EnumSet.of(ModeOfInheritance.AUTOSOMAL_DOMINANT));
        ve.setContributesToGeneScoreUnderMode(ModeOfInheritance.AUTOSOMAL_DOMINANT);
        if (predictor != null)
            ve.setPathogenicityData(PathogenicityData.of(PathogenicityScore.of(PathogenicitySource.REVEL, predictor)));
        return ve;
    }
}
