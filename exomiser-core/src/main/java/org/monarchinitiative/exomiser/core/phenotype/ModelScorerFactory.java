package org.monarchinitiative.exomiser.core.phenotype;

import java.util.List;

/**
 * Factory for creating {@link ModelScorer} instances given the phenotype query context.
 *
 * <p>Two factory methods correspond to the two scoring modes used by Exomiser prioritisers:
 * <ul>
 *   <li>{@link #forSingleCrossSpecies} — single-species cross comparison (e.g. HP vs MP in Phive)</li>
 *   <li>{@link #forMultiCrossSpecies} — multi-species comparison normalised against a human reference (HiPhive)</li>
 * </ul>
 *
 * <p>Classical implementations (Phenodigm) use the {@link PhenotypeMatcher} and
 * {@link QueryPhenotypeMatch} arguments for IC-based scoring.  Neural implementations
 * (e.g. INDIGENA) typically ignore those arguments and operate purely on embedding
 * vectors derived from {@code queryTerms}.
 *
 * <p>To add a new phenotype similarity method, implement this interface and wire it
 * via Spring configuration — no changes to the prioritisers themselves are required.
 */
public interface ModelScorerFactory {

    /**
     * Create a scorer for a single-species cross comparison (e.g. HP→MP for mouse models).
     * Used by {@code PhivePriority}.
     *
     * @param queryTerms         patient phenotype terms
     * @param phenotypeMatcher   pre-computed best-match data for the target organism
     */
    <T extends Model> ModelScorer<T> forSingleCrossSpecies(
            List<PhenotypeTerm> queryTerms,
            PhenotypeMatcher phenotypeMatcher);

    /**
     * Create a scorer for multi-species comparison, normalised against the human HP-HP reference.
     * Used by {@code HiPhivePriority}.
     *
     * @param queryTerms                     patient phenotype terms
     * @param referenceQueryPhenotypeMatch   HP-HP self-match used for score normalisation
     * @param phenotypeMatcher               pre-computed best-match data for the target organism
     */
    <T extends Model> ModelScorer<T> forMultiCrossSpecies(
            List<PhenotypeTerm> queryTerms,
            QueryPhenotypeMatch referenceQueryPhenotypeMatch,
            PhenotypeMatcher phenotypeMatcher);
}
