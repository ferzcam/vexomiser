package org.monarchinitiative.exomiser.core.phenotype;

import org.monarchinitiative.exomiser.core.prioritisers.model.GeneModel;
import org.monarchinitiative.exomiser.core.model.Gene;
import java.util.List;
import java.util.Map;
import java.util.HashMap;
import java.util.Set;
import java.util.ArrayList;
import java.util.Collections;
import java.util.stream.Collectors;

/** Per-analysis native adapter for a precomputed EmbedPVP case entity. */
public final class TransdModelScorerFactory implements LearnedModelScorerFactory {
    private final TransdCheckpoint checkpoint;
    private final String caseId;
    private final Map<Integer, Double> normalized;

    public TransdModelScorerFactory(TransdCheckpoint checkpoint, String caseId) {
        this.checkpoint = checkpoint;
        if (caseId == null || caseId.isBlank()) throw new IllegalArgumentException("TransD caseId is required");
        this.caseId = caseId;
        this.normalized = null;
    }

    private TransdModelScorerFactory(TransdCheckpoint checkpoint, String caseId, Map<Integer, Double> normalized) {
        this.checkpoint = checkpoint;
        this.caseId = caseId;
        this.normalized = Map.copyOf(normalized);
    }

    @Override
    public ModelScorerFactory withCandidateGenes(List<Gene> genes) {
        Set<Integer> ids = genes.stream().map(Gene::entrezGeneId).collect(Collectors.toSet());
        if (!checkpoint.hasCase(caseId)) throw new IllegalArgumentException("TransD bundle has no case entity: " + caseId);
        Map<Integer, Double> raw = new HashMap<>();
        for (int id : ids) {
            if (checkpoint.hasGene(Integer.toString(id)))
                raw.put(id, checkpoint.scoreGeneCase(Integer.toString(id), caseId));
        }
        int missingCount = ids.size() - raw.size();
        org.slf4j.LoggerFactory.getLogger(TransdModelScorerFactory.class).info(
                "TransD native case {}: {} covered genes, {} absent genes in analysis context", caseId, raw.size(), missingCount);
        double median = median(new ArrayList<>(raw.values()));
        for (int id : ids) raw.putIfAbsent(id, median);
        if (raw.isEmpty()) return new TransdModelScorerFactory(checkpoint, caseId, Map.of());
        double min = raw.values().stream().mapToDouble(Double::doubleValue).min().orElseThrow();
        double max = raw.values().stream().mapToDouble(Double::doubleValue).max().orElseThrow();
        Map<Integer, Double> scaled = new HashMap<>();
        raw.forEach((id, score) -> scaled.put(id, max == min ? 0.0 : (score - min) / (max - min)));
        return new TransdModelScorerFactory(checkpoint, caseId, scaled);
    }

    private static double median(List<Double> values) {
        if (values.isEmpty()) return 0.0;
        Collections.sort(values);
        int n = values.size();
        return n % 2 == 1 ? values.get(n / 2) : (values.get(n / 2 - 1) + values.get(n / 2)) / 2;
    }

    @Override
    public <T extends Model> ModelScorer<T> forSingleCrossSpecies(List<PhenotypeTerm> queryTerms, PhenotypeMatcher matcher) {
        return scorer();
    }

    @Override
    public <T extends Model> ModelScorer<T> forMultiCrossSpecies(List<PhenotypeTerm> queryTerms,
            QueryPhenotypeMatch reference, PhenotypeMatcher matcher) {
        return scorer();
    }

    private <T extends Model> ModelScorer<T> scorer() {
        if (normalized == null) throw new IllegalStateException("TransD scorer requires candidate gene context");
        return model -> {
            if (!(model instanceof GeneModel geneModel))
                throw new IllegalArgumentException("TransD native phenotype scoring requires GeneModel");
            Double compatibility = normalized.get(geneModel.entrezGeneId());
            if (compatibility == null) throw new IllegalArgumentException("TransD model gene is outside analysis context: " + geneModel.entrezGeneId());
            return ModelPhenotypeMatch.of(compatibility, model, List.of());
        };
    }
}
