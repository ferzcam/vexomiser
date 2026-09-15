package org.monarchinitiative.exomiser.core.analysis;

import org.monarchinitiative.exomiser.api.v1.AnalysisProto;
import org.monarchinitiative.exomiser.core.analysis.score.WeightedVariantScorer;

import java.nio.file.Path;

/** Per-analysis, validated configuration for native phenotype and independent variant scoring. */
public record LearnedScoringOptions(Method method, Path embeddings, Path transdBundle, String caseId,
                                    Weighted weighted) {
    public enum Method { PHENODIGM, INDIGENA, EMBEDPVP_TRANSD }
    public record Weighted(double phenotypeWeight, Path genePhenotypes, Path patientHpo,
                           WeightedVariantScorer.MissingPolicy pathogenicityMissing,
                           WeightedVariantScorer.MissingPolicy phenotypeMissing) {}

    public static LearnedScoringOptions defaults() {
        return new LearnedScoringOptions(Method.PHENODIGM, null, null, "", null);
    }

    public static LearnedScoringOptions from(AnalysisProto.Analysis proto) {
        var scorer = proto.getPhenotypeScorer();
        String methodText = scorer.getMethod().isBlank() ? "PHENODIGM" : scorer.getMethod();
        Method method;
        try { method = Method.valueOf(methodText); }
        catch (IllegalArgumentException e) { throw new IllegalArgumentException("Unknown phenotype scorer: " + methodText, e); }
        Path embeddings = path(scorer.getEmbeddings());
        Path bundle = path(scorer.getTransdBundle());
        String caseId = scorer.getCaseId();
        switch (method) {
            case PHENODIGM -> {
                if (embeddings != null || bundle != null || !caseId.isBlank())
                    throw new IllegalArgumentException("PHENODIGM does not accept learned model resources");
            }
            case INDIGENA -> {
                if (embeddings == null || bundle != null)
                    throw new IllegalArgumentException("INDIGENA requires embeddings and no TransD bundle");
            }
            case EMBEDPVP_TRANSD -> {
                if (bundle == null || caseId.isBlank() || embeddings != null)
                    throw new IllegalArgumentException("EMBEDPVP_TRANSD requires transdBundle and caseId only");
            }
        }
        Weighted weighted = null;
        if (proto.hasWeightedVariants()) {
            var options = proto.getWeightedVariants();
            if (method == Method.PHENODIGM)
                throw new IllegalArgumentException("weightedVariants requires INDIGENA or EMBEDPVP_TRANSD");
            Path annotations = path(options.getGenePhenotypes());
            if (method == Method.INDIGENA && annotations == null)
                throw new IllegalArgumentException("INDIGENA weightedVariants requires genePhenotypes");
            if (caseId.isBlank())
                throw new IllegalArgumentException("weightedVariants requires phenotypeScorer.caseId for output provenance");
            if (!options.hasPhenotypeWeight())
                throw new IllegalArgumentException("weightedVariants requires an explicit phenotypeWeight (zero is valid)");
            double weight = options.getPhenotypeWeight();
            if (!Double.isFinite(weight) || weight < 0 || weight > 1)
                throw new IllegalArgumentException("phenotypeWeight must be within [0,1]");
            weighted = new Weighted(weight, annotations, path(options.getPatientHpo()),
                    missing(options.getPathogenicityMissing(), WeightedVariantScorer.MissingPolicy.ZERO),
                    missing(options.getPhenotypeMissing(), WeightedVariantScorer.MissingPolicy.MEDIAN));
        }
        return new LearnedScoringOptions(method, embeddings, bundle, caseId, weighted);
    }

    private static Path path(String text) { return text.isBlank() ? null : Path.of(text); }
    private static WeightedVariantScorer.MissingPolicy missing(String text, WeightedVariantScorer.MissingPolicy fallback) {
        if (text.isBlank()) return fallback;
        try { return WeightedVariantScorer.MissingPolicy.valueOf(text); }
        catch (IllegalArgumentException e) { throw new IllegalArgumentException("Unknown missing policy: " + text, e); }
    }
}
