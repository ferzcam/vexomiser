package org.monarchinitiative.exomiser.core.phenotype;

import java.util.ArrayList;
import java.util.List;

/**
 * {@link ModelScorer} that computes phenotype similarity using INDIGENA embedding vectors
 * instead of the classical Phenodigm IC-weighted semantic similarity.
 *
 * <p>Similarity between two phenotype sets is the Best-Match Average (BMA) of
 * sigmoid(dot-product) similarities between their embedding vectors, averaged
 * symmetrically (disease-centric and gene-centric directions).
 *
 * <p>Terms absent from the embedding vocabulary are silently skipped; if neither
 * the query nor the model has any known terms the score is 0.
 */
public class IndigenaModelScorer<T extends Model> implements ModelScorer<T> {

    private final IndigenaEmbeddings embeddings;
    // Pre-looked-up query embedding vectors (one per matched query phenotype).
    private final float[][] queryVecs;

    /**
     * @param queryPhenotypeTerms patient phenotype terms (HP IDs like {@code HP:0001234})
     * @param embeddings          loaded INDIGENA embedding store
     */
    public IndigenaModelScorer(List<PhenotypeTerm> queryPhenotypeTerms, IndigenaEmbeddings embeddings) {
        this.embeddings = embeddings;
        List<float[]> vecs = new ArrayList<>(queryPhenotypeTerms.size());
        for (PhenotypeTerm term : queryPhenotypeTerms) {
            float[] v = embeddings.getByPhenotypeId(term.id());
            if (v != null) vecs.add(v);
        }
        this.queryVecs = vecs.toArray(new float[0][]);
    }

    @Override
    public ModelPhenotypeMatch<T> scoreModel(T model) {
        if (queryVecs.length == 0) {
            return ModelPhenotypeMatch.of(0.0, model, List.of());
        }

        List<float[]> modelVecs = new ArrayList<>();
        for (String phenoId : model.phenotypeIds()) {
            float[] v = embeddings.getByPhenotypeId(phenoId);
            if (v != null) modelVecs.add(v);
        }

        if (modelVecs.isEmpty()) {
            return ModelPhenotypeMatch.of(0.0, model, List.of());
        }

        float[][] mv = modelVecs.toArray(new float[0][]);
        double score = bma(queryVecs, mv);
        return ModelPhenotypeMatch.of(score, model, List.of());
    }

    /**
     * Symmetric BMA: average of disease-centric and gene-centric best-match averages.
     *
     * @param q query (patient) phenotype vectors  (Q x D)
     * @param m model (gene/disease) phenotype vectors (M x D)
     */
    private double bma(float[][] q, float[][] m) {
        int Q = q.length;
        int M = m.length;

        // Disease-centric: for each query term, best match across model terms.
        double diseaseCentric = 0;
        for (float[] qVec : q) {
            float best = 0;
            for (float[] mVec : m) best = Math.max(best, sigmoid(dot(qVec, mVec)));
            diseaseCentric += best;
        }
        diseaseCentric /= Q;

        // Gene-centric: for each model term, best match across query terms.
        double geneCentric = 0;
        for (float[] mVec : m) {
            float best = 0;
            for (float[] qVec : q) best = Math.max(best, sigmoid(dot(qVec, mVec)));
            geneCentric += best;
        }
        geneCentric /= M;

        return (diseaseCentric + geneCentric) / 2.0;
    }

    private float dot(float[] a, float[] b) {
        float sum = 0;
        for (int i = 0; i < a.length; i++) sum += a[i] * b[i];
        return sum;
    }

    private float sigmoid(float x) {
        return 1.0f / (1.0f + (float) Math.exp(-x));
    }
}
