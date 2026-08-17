package org.monarchinitiative.exomiser.core.phenotype;

import org.monarchinitiative.exomiser.core.prioritisers.model.GeneModel;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * {@link ModelScorer} implementing EmbedPVP's phenotype term (its "MS").
 *
 * <p>The authoritative definition is
 * {@code EmbedPVP/embedpvp/main.py::calculate_graph_embed_similarity}:
 *
 * <pre>
 *     MS(gene) = (cosine(e_case, e_gene) + 1) / 2
 * </pre>
 *
 * <p>and the variant's score is the max of that over the variant's genes. Note what this
 * is <em>not</em>: unlike {@link IndigenaModelScorer}, no phenotype term sets are compared.
 * The model's {@link Model#phenotypeIds()} are ignored entirely; the only thing read off
 * the model is its Entrez gene id. Every model of a gene therefore receives the same score,
 * and the prioritiser's max-over-models reduces to MS(gene) — which is what EmbedPVP means.
 *
 * <p>Because the score depends on the case rather than on the query terms, the case
 * identity must be supplied out of band; see {@link EmbedpvpModelScorerFactory} and
 * {@code docs/embedpvp_scorer.md}.
 *
 * <p><b>Known gap.</b> The prioritiser only calls a scorer for genes that have at least
 * one model. EmbedPVP's MS is defined for every gene with an embedding, so genes with no
 * mouse ortholog and no disease association score 0 here where EmbedPVP would give them a
 * real value. Reproducing EmbedPVP over that slice needs a {@code Prioritiser}, not a
 * {@code ModelScorer}; the variant-level harness takes that route instead.
 */
public class EmbedpvpModelScorer<T extends Model> implements ModelScorer<T> {

    private static final Logger logger = LoggerFactory.getLogger(EmbedpvpModelScorer.class);

    private final EmbedpvpEmbeddings embeddings;
    private final String caseId;
    private final float[] caseVector;

    /**
     * @param caseId     the case being scored, e.g. {@code PAVS:A0000003}
     * @param queryTerms the patient's phenotype terms — used only to check that
     *                   {@code caseId} really is this patient
     * @param embeddings loaded gene and case embeddings
     */
    public EmbedpvpModelScorer(String caseId, List<PhenotypeTerm> queryTerms,
                               EmbedpvpEmbeddings embeddings) {
        this.embeddings = embeddings;
        this.caseId = caseId;
        this.caseVector = embeddings.getCaseVector(caseId);
        if (caseVector == null) {
            logger.warn("No EmbedPVP case embedding for '{}' — every gene will score 0. "
                    + "Check the case-embedding export covers this case.", caseId);
        } else {
            warnIfPhenotypesDisagree(queryTerms);
        }
    }

    /**
     * Guard against a case embedding attached to the wrong patient. An exact match is not
     * expected — Exomiser drops query terms absent from its own HPO release — so only a
     * completely empty intersection is reported.
     */
    private void warnIfPhenotypesDisagree(List<PhenotypeTerm> queryTerms) {
        Set<String> exported = embeddings.getCasePhenotypes(caseId);
        if (exported.isEmpty() || queryTerms == null || queryTerms.isEmpty()) {
            return;
        }
        Set<String> query = new HashSet<>();
        for (PhenotypeTerm term : queryTerms) {
            query.add(term.id());
        }
        query.retainAll(exported);
        if (query.isEmpty()) {
            logger.warn("EmbedPVP case '{}' was exported with phenotypes {} but Exomiser's "
                            + "query terms share none of them — the case embedding may be "
                            + "wired to the wrong patient.", caseId, exported);
        }
    }

    @Override
    public ModelPhenotypeMatch<T> scoreModel(T model) {
        if (caseVector == null || !(model instanceof GeneModel geneModel)) {
            return ModelPhenotypeMatch.of(0.0, model, List.of());
        }
        float[] geneVector = embeddings.getGeneVector(geneModel.entrezGeneId());
        if (geneVector == null) {
            // EmbedPVP's own fallback: a gene outside the embedding vocabulary contributes
            // nothing, i.e. max_gene_score = 0.
            return ModelPhenotypeMatch.of(0.0, model, List.of());
        }
        return ModelPhenotypeMatch.of(similarity(caseVector, geneVector), model, List.of());
    }

    /** {@code (cosine + 1) / 2}, so the score lands on [0, 1] as EmbedPVP expects. */
    static double similarity(float[] a, float[] b) {
        double dot = 0;
        double normA = 0;
        double normB = 0;
        int n = Math.min(a.length, b.length);
        for (int i = 0; i < n; i++) {
            dot += (double) a[i] * b[i];
            normA += (double) a[i] * a[i];
            normB += (double) b[i] * b[i];
        }
        if (normA == 0 || normB == 0) {
            return 0.0;
        }
        return (dot / (Math.sqrt(normA) * Math.sqrt(normB)) + 1.0) / 2.0;
    }
}
