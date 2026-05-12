package org.monarchinitiative.exomiser.core.phenotype;

import java.util.List;

/**
 * {@link ModelScorerFactory} that scores gene models using pure Resnik IC BMA
 * instead of the Phenodigm {@code sqrt(IC × simJ)} formula.
 *
 * <p>The same HP-HP and HP-MP database mappings used by HiPhive and Phive are
 * reused — only the pairwise similarity function changes.  This makes
 * Resnik-HiPhive / Resnik-Phive / Resnik-PhenIX a controlled comparison against
 * INDIGENA-HiPhive etc., isolating the effect of the embedding space from all
 * other pipeline differences.
 */
public class ResnikModelScorerFactory implements ModelScorerFactory {

    public ResnikModelScorerFactory() {}

    @Override
    public <T extends Model> ModelScorer<T> forSingleCrossSpecies(
            List<PhenotypeTerm> queryTerms,
            PhenotypeMatcher phenotypeMatcher) {
        return new ResnikModelScorer<>(phenotypeMatcher);
    }

    @Override
    public <T extends Model> ModelScorer<T> forMultiCrossSpecies(
            List<PhenotypeTerm> queryTerms,
            QueryPhenotypeMatch referenceQueryPhenotypeMatch,
            PhenotypeMatcher phenotypeMatcher) {
        return new ResnikModelScorer<>(phenotypeMatcher);
    }
}
