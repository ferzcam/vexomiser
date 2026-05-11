package org.monarchinitiative.exomiser.core.phenotype;

import java.util.List;

/**
 * {@link ModelScorerFactory} that delegates to the classical Phenodigm IC-based
 * semantic similarity algorithm ({@link PhenodigmModelScorer}).
 *
 * <p>This is the default factory used by Exomiser's existing prioritisers.
 * The {@code queryTerms} parameter is not used directly here because Phenodigm
 * pre-computes all pairwise term similarities into the {@link PhenotypeMatcher}.
 */
public class PhenodigmModelScorerFactory implements ModelScorerFactory {

    public static final PhenodigmModelScorerFactory INSTANCE = new PhenodigmModelScorerFactory();

    private PhenodigmModelScorerFactory() {}

    @Override
    public <T extends Model> ModelScorer<T> forSingleCrossSpecies(
            List<PhenotypeTerm> queryTerms,
            PhenotypeMatcher phenotypeMatcher) {
        return PhenodigmModelScorer.forSingleCrossSpecies(phenotypeMatcher);
    }

    @Override
    public <T extends Model> ModelScorer<T> forMultiCrossSpecies(
            List<PhenotypeTerm> queryTerms,
            QueryPhenotypeMatch referenceQueryPhenotypeMatch,
            PhenotypeMatcher phenotypeMatcher) {
        return PhenodigmModelScorer.forMultiCrossSpecies(referenceQueryPhenotypeMatch, phenotypeMatcher);
    }
}
