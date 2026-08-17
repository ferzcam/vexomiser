package org.monarchinitiative.exomiser.core.phenotype;

import java.util.List;

/**
 * {@link ModelScorerFactory} producing {@link EmbedpvpModelScorer}s.
 *
 * <p>EmbedPVP's phenotype score is a similarity between a <em>case</em> embedding and a
 * <em>gene</em> embedding, so the scorer needs to know which case it is scoring — something
 * the {@link ModelScorerFactory} interface does not pass. The case id therefore arrives out
 * of band, in this precedence order:
 *
 * <ol>
 *   <li>the {@code caseId} constructor argument (JPype harnesses, which score many cases in
 *       one JVM, use {@link #withCase(String)} to rebind it);</li>
 *   <li>the {@code embedpvp.caseId} system property, when the constructor argument is null
 *       (the CLI path, which forks one JVM per case).</li>
 * </ol>
 *
 * <p>Both factory methods ignore their {@link PhenotypeMatcher} and
 * {@link QueryPhenotypeMatch} arguments: EmbedPVP does no cross-species term matching, and
 * no per-organism normalisation applies to an entity-level cosine.
 *
 * <p>See {@code docs/embedpvp_scorer.md} for the export contract.
 */
public class EmbedpvpModelScorerFactory implements ModelScorerFactory {

    /** System property consulted when no case id is supplied explicitly. */
    public static final String CASE_ID_PROPERTY = "embedpvp.caseId";

    private final EmbedpvpEmbeddings embeddings;
    private final String caseId;

    public EmbedpvpModelScorerFactory(EmbedpvpEmbeddings embeddings) {
        this(embeddings, null);
    }

    public EmbedpvpModelScorerFactory(EmbedpvpEmbeddings embeddings, String caseId) {
        this.embeddings = embeddings;
        this.caseId = caseId != null ? caseId : System.getProperty(CASE_ID_PROPERTY);
    }

    /** A factory over the same embeddings, bound to a different case. */
    public EmbedpvpModelScorerFactory withCase(String otherCaseId) {
        return new EmbedpvpModelScorerFactory(embeddings, otherCaseId);
    }

    public String caseId() {
        return caseId;
    }

    @Override
    public <T extends Model> ModelScorer<T> forSingleCrossSpecies(
            List<PhenotypeTerm> queryTerms,
            PhenotypeMatcher phenotypeMatcher) {
        return new EmbedpvpModelScorer<>(caseId, queryTerms, embeddings);
    }

    @Override
    public <T extends Model> ModelScorer<T> forMultiCrossSpecies(
            List<PhenotypeTerm> queryTerms,
            QueryPhenotypeMatch referenceQueryPhenotypeMatch,
            PhenotypeMatcher phenotypeMatcher) {
        return new EmbedpvpModelScorer<>(caseId, queryTerms, embeddings);
    }
}
