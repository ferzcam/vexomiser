package org.monarchinitiative.exomiser.core.phenotype;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.math.BigDecimal;
import java.math.RoundingMode;

/** INDIGENA BMA over benchmark pooled gene-to-MP annotations, independent of database models. */
public final class PooledGenePhenotypeScorer {
    private final IndigenaEmbeddings embeddings;
    private final Map<String, List<String>> annotations;

    public PooledGenePhenotypeScorer(IndigenaEmbeddings embeddings, Path csv) throws IOException {
        this.embeddings = embeddings;
        Map<String, List<String>> parsed = new HashMap<>();
        try (var reader = Files.newBufferedReader(csv)) {
            String line = reader.readLine(); // header, as benchmark pairs() does
            if (line == null) throw new IOException("Empty gene phenotype CSV: " + csv);
            int number = 1;
            while ((line = reader.readLine()) != null) {
                number++;
                String[] fields = line.split(",", -1);
                if (fields.length != 2 || fields[0].isBlank() || fields[1].isBlank())
                    throw new IOException(csv + ":" + number + " requires gene,phenotype pair");
                parsed.computeIfAbsent(geneId(fields[0]), ignored -> new ArrayList<>()).add(fields[1]);
            }
        }
        this.annotations = Map.copyOf(parsed);
    }

    public Map<String, Double> score(Set<String> geneIds, List<String> patientTerms) {
        List<float[]> query = vectors(patientTerms);
        Map<String, Double> scores = new HashMap<>();
        for (String geneId : geneIds) {
            List<float[]> gene = vectors(annotations.getOrDefault(geneId, List.of()));
            double score = query.isEmpty() || gene.isEmpty() ? 0.0 : bma(query, gene);
            // The reference phenotype table formats each score to six decimals before combination.
            scores.put(geneId, sixDecimal(score));
        }
        return scores;
    }

    static double sixDecimal(double score) {
        if (!Double.isFinite(score)) throw new IllegalArgumentException("Non-finite pooled phenotype score");
        return new BigDecimal(score).setScale(6, RoundingMode.HALF_EVEN).doubleValue();
    }

    private List<float[]> vectors(List<String> ids) {
        List<float[]> vectors = new ArrayList<>();
        for (String id : ids) {
            float[] vector = id.startsWith("http://") || id.startsWith("https://")
                    ? embeddings.getByIri(id) : embeddings.getByPhenotypeId(id);
            if (vector != null) vectors.add(vector);
        }
        return vectors;
    }

    private static String geneId(String text) {
        int slash = text.lastIndexOf('/');
        return slash < 0 ? text : text.substring(slash + 1);
    }

    private static double bma(List<float[]> query, List<float[]> gene) {
        double queryCentric = 0;
        for (float[] q : query) {
            double best = 0;
            for (float[] g : gene) best = Math.max(best, similarity(q, g));
            queryCentric += best;
        }
        double geneCentric = 0;
        for (float[] g : gene) {
            double best = 0;
            for (float[] q : query) best = Math.max(best, similarity(q, g));
            geneCentric += best;
        }
        return (queryCentric / query.size() + geneCentric / gene.size()) / 2;
    }

    private static double similarity(float[] a, float[] b) {
        float dot = 0;
        for (int i = 0; i < a.length; i++) dot += a[i] * b[i];
        return 1.0 / (1.0 + Math.exp(-dot));
    }
}
