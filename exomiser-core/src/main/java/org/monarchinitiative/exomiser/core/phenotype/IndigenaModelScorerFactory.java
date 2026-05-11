package org.monarchinitiative.exomiser.core.phenotype;

import java.util.List;

/**
 * {@link ModelScorerFactory} that uses INDIGENA embedding-based BMA similarity
 * instead of Phenodigm IC-weighted semantic similarity.
 *
 * <p>Both factory methods return an {@link IndigenaModelScorer} driven solely by
 * {@code queryTerms} and the pre-loaded {@link IndigenaEmbeddings}.  The
 * {@link PhenotypeMatcher} and {@link QueryPhenotypeMatch} arguments are ignored
 * because the UPheno embedding space already covers HP, MP, and related ontologies
 * across species — no separate cross-species mapping step is needed.
 *
 * <p>Usage: load embeddings once at application startup and inject this factory
 * wherever a {@link ModelScorerFactory} is expected.
 */
public class IndigenaModelScorerFactory implements ModelScorerFactory {

    private final IndigenaEmbeddings embeddings;

    public IndigenaModelScorerFactory(IndigenaEmbeddings embeddings) {
        this.embeddings = embeddings;
    }

    @Override
    public <T extends Model> ModelScorer<T> forSingleCrossSpecies(
            List<PhenotypeTerm> queryTerms,
            PhenotypeMatcher phenotypeMatcher) {
        return new IndigenaModelScorer<>(queryTerms, embeddings);
    }

    @Override
    public <T extends Model> ModelScorer<T> forMultiCrossSpecies(
            List<PhenotypeTerm> queryTerms,
            QueryPhenotypeMatch referenceQueryPhenotypeMatch,
            PhenotypeMatcher phenotypeMatcher) {
        return new IndigenaModelScorer<>(queryTerms, embeddings);
    }
}
