package org.monarchinitiative.exomiser.core.phenotype;

import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * {@link ModelScorer} that computes phenotype similarity using pure Resnik IC
 * as the pairwise measure with symmetric BMA as the groupwise aggregation.
 *
 * <p>Pairwise similarity between two terms is {@code IC(LCS(t1, t2))} — the
 * information content of the lowest common subsumer — taken directly from the
 * {@code ic} field of the pre-computed {@link PhenotypeMatch} entries provided
 * by the {@link PhenotypeMatcher}.  The same HP-HP / HP-MP database mappings
 * used by Phenodigm (HiPhive / Phive) are reused here; only the scoring
 * formula changes from {@code sqrt(IC × simJ)} to {@code IC}.
 *
 * <p>Groupwise score: symmetric BMA — average of:
 * <ul>
 *   <li>query-centric: mean over query terms of best IC match in model</li>
 *   <li>model-centric: mean over model terms of best IC match in query</li>
 * </ul>
 *
 * <p>Scores are not normalised against a theoretical maximum; only relative
 * rankings matter for evaluation.
 */
public class ResnikModelScorer<T extends Model> implements ModelScorer<T> {

    private final PhenotypeMatcher phenotypeMatcher;

    ResnikModelScorer(PhenotypeMatcher phenotypeMatcher) {
        this.phenotypeMatcher = phenotypeMatcher;
    }

    @Override
    public ModelPhenotypeMatch<T> scoreModel(T model) {
        List<String> modelTermIds = model.phenotypeIds();
        if (modelTermIds.isEmpty()) {
            return ModelPhenotypeMatch.of(0.0, model, List.of());
        }

        Map<PhenotypeTerm, Set<PhenotypeMatch>> termMatches = phenotypeMatcher.getTermPhenotypeMatches();
        if (termMatches.isEmpty()) {
            return ModelPhenotypeMatch.of(0.0, model, List.of());
        }

        Set<String> modelTermSet = new HashSet<>(modelTermIds);

        // Query-centric: for each query term, find best IC match among model terms.
        // Also accumulate model-centric best ICs as a by-product.
        double queryCentricSum = 0.0;
        Map<String, Double> modelTermBestIC = new HashMap<>();

        for (Map.Entry<PhenotypeTerm, Set<PhenotypeMatch>> entry : termMatches.entrySet()) {
            double bestIC = 0.0;
            for (PhenotypeMatch match : entry.getValue()) {
                if (modelTermSet.contains(match.matchPhenotypeId())) {
                    double ic = match.ic();
                    if (ic > bestIC) bestIC = ic;
                    modelTermBestIC.merge(match.matchPhenotypeId(), ic, Math::max);
                }
            }
            queryCentricSum += bestIC;
        }

        double queryCentric = queryCentricSum / termMatches.size();

        // Model-centric: mean over model terms of best IC match from any query term.
        double modelCentricSum = 0.0;
        for (String modelTermId : modelTermIds) {
            modelCentricSum += modelTermBestIC.getOrDefault(modelTermId, 0.0);
        }
        double modelCentric = modelCentricSum / modelTermIds.size();

        double score = (queryCentric + modelCentric) / 2.0;
        return ModelPhenotypeMatch.of(score, model, List.of());
    }
}
