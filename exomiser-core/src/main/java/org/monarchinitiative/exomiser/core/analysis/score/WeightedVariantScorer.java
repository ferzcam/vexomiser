/*
 * The Exomiser - A tool to annotate and prioritize genomic variants
 *
 * Copyright (c) 2016-2021 Queen Mary University of London.
 * Copyright (c) 2012-2016 Charité Universitätsmedizin Berlin and Genome Research Ltd.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package org.monarchinitiative.exomiser.core.analysis.score;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Set;

/**
 * Pure per-case variant scoring engine that combines a per-gene phenotype score (MS) with a
 * per-variant pathogenicity score (GP) into the EmbedPVP combined score
 *
 * <pre>S(variant) = w * MS(variant) + (1 - w) * GP(variant)</pre>
 *
 * where MS(variant) is the maximum phenotype score over the variant's genes and both component
 * columns are min-max normalised within the case before combining. This mirrors the
 * variant-level combination semantics of the EmbedPVP2 evaluation harness
 * ({@code code/eval/combine_variant_scores.py}) for the Table 3 setup: CADD pathogenicity with
 * the ZERO missing policy and phenotype scores with the MEDIAN missing policy.
 *
 * The engine is stateless and deterministic. Input rows are merged by genomic key (duplicate
 * rows union their gene sets and keep the maximum non-missing pathogenicity score); an
 * optional {@code chr} prefix on the contig is canonicalised the same way as in the reference
 * harness. Alleles and positions are used as given and are not biologically normalised.
 * Non-finite scores or weights are rejected rather than silently ranked.
 */
public final class WeightedVariantScorer {

    /**
     * Policy for imputing a missing component value before per-case min-max normalisation.
     */
    public enum MissingPolicy {
        /** Impute 0.0 (the EmbedPVP {@code fillna(0)} behaviour). */
        ZERO,
        /** Impute the median over the distinct variants of the case that do have a value. */
        MEDIAN
    }

    /**
     * Genomic identity of a variant. The optional {@code chr} prefix on the contig is
     * canonicalised away so that, e.g., {@code chr1} and {@code 1} denote the same contig.
     */
    public record VariantKey(String contig, int start, String ref, String alt) implements Comparable<VariantKey> {

        public VariantKey {
            Objects.requireNonNull(contig, "contig must not be null");
            Objects.requireNonNull(ref, "ref must not be null");
            Objects.requireNonNull(alt, "alt must not be null");
            if (contig.startsWith("chr")) {
                contig = contig.substring(3);
            }
        }

        @Override
        public int compareTo(VariantKey other) {
            int c = contig.compareTo(other.contig);
            if (c != 0) {
                return c;
            }
            c = Integer.compare(start, other.start);
            if (c != 0) {
                return c;
            }
            c = ref.compareTo(other.ref);
            if (c != 0) {
                return c;
            }
            return alt.compareTo(other.alt);
        }
    }

    /**
     * One input row: a variant with an optional pathogenicity score and its gene assignments.
     *
     * <p>{@code pathogenicityScore} may be {@code null} for a variant the predictor could not
     * score (missing is distinct from an explicit zero). It must be finite when present.
     * Duplicate rows (same genomic key) are merged by the scorer.</p>
     */
    public record VariantInput(VariantKey key, Double pathogenicityScore, Set<String> geneIds) {

        public VariantInput {
            Objects.requireNonNull(key, "key must not be null");
            Objects.requireNonNull(geneIds, "geneIds must not be null");
            geneIds = Set.copyOf(geneIds);
            if (pathogenicityScore != null && !Double.isFinite(pathogenicityScore)) {
                throw new IllegalArgumentException(
                        "pathogenicityScore must be finite or null (missing), got " + pathogenicityScore);
            }
        }

        /**
         * @return true if this row carries no pathogenicity score.
         */
        public boolean pathogenicityMissing() {
            return pathogenicityScore == null;
        }
    }

    /**
     * Scoring parameters: phenotype weight in the combined score and the imputation policies
     * for missing component values.
     */
    public record Settings(double phenotypeWeight, MissingPolicy pathogenicityMissing, MissingPolicy phenotypeMissing) {

        public Settings {
            Objects.requireNonNull(pathogenicityMissing, "pathogenicityMissing must not be null");
            Objects.requireNonNull(phenotypeMissing, "phenotypeMissing must not be null");
            if (!Double.isFinite(phenotypeWeight) || phenotypeWeight < 0.0 || phenotypeWeight > 1.0) {
                throw new IllegalArgumentException(
                        "phenotypeWeight must be finite and within [0, 1], got " + phenotypeWeight);
            }
        }

        /**
         * Table 3 defaults: CADD pathogenicity missing values imputed as zero, phenotype
         * missing values imputed as the median over scoreable variants of the case.
         */
        public static Settings table3(double phenotypeWeight) {
            return new Settings(phenotypeWeight, MissingPolicy.ZERO, MissingPolicy.MEDIAN);
        }
    }

    /**
     * Score for one distinct variant of a case.
     *
     * @param key                         genomic key (contig canonicalised)
     * @param geneIds                     union of gene ids across all input rows for this variant
     * @param rawPhenotypeScore           raw phenotype score before imputation (max over the
     *                                    variant's genes present in the score table);
     *                                    {@code null} when no gene of the variant was scoreable
     * @param rawPathogenicityScore       raw pathogenicity score before imputation (max of the
     *                                    non-missing values across duplicate rows);
     *                                    {@code null} when all rows were missing
     * @param phenotypeMissing            true when no phenotype score was available for this variant
     * @param pathogenicityMissing        true when no pathogenicity score was available for this variant
     * @param imputedPhenotypeScore       raw phenotype score after applying the missing policy
     * @param imputedPathogenicityScore   raw pathogenicity score after applying the missing policy
     * @param normalizedPhenotypeScore    min-max normalised phenotype score in [0, 1] (0 for a constant column)
     * @param normalizedPathogenicityScore min-max normalised pathogenicity score in [0, 1] (0 for a constant column)
     * @param combinedScore               {@code w * normalizedPhenotypeScore + (1 - w) * normalizedPathogenicityScore}
     * @param rank                        average tie rank: number of strictly greater scores + (number of exact ties + 1) / 2
     * @param rankStrict                   number of strictly greater scores + 1
     */
    public record VariantScore(
            VariantKey key,
            Set<String> geneIds,
            Double rawPhenotypeScore,
            Double rawPathogenicityScore,
            boolean phenotypeMissing,
            boolean pathogenicityMissing,
            double imputedPhenotypeScore,
            double imputedPathogenicityScore,
            double normalizedPhenotypeScore,
            double normalizedPathogenicityScore,
            double combinedScore,
            double rank,
            double rankStrict) {

        public VariantScore {
            Objects.requireNonNull(key, "key must not be null");
            Objects.requireNonNull(geneIds, "geneIds must not be null");
            geneIds = Set.copyOf(geneIds);
        }
    }

    /**
     * Score all variants of one case and return them sorted by descending combined score.
     *
     * <p>Semantics (matching the reference harness):</p>
     * <ul>
     * <li>input rows are merged by genomic key: gene ids are unioned and the maximum
     * non-missing pathogenicity score is kept (missing is distinct from an explicit zero);</li>
     * <li>the phenotype score of a variant is the maximum over the scores of its genes present
     * in {@code genePhenotypeScores}; an explicit zero score is a real score; a variant with no
     * scoreable gene is missing;</li>
     * <li>missing values are imputed per {@link Settings}: ZERO imputes 0.0, MEDIAN imputes the
     * median over the distinct variants of the case that are not missing (0.0 when none are);</li>
     * <li>both component columns are min-max normalised per case over the imputed values; a
     * constant column maps to 0;</li>
     * <li>the combined score is {@code w * MS + (1 - w) * GP}; exact score equality defines
     * ties, which receive the average tie rank.</li>
     * </ul>
     *
     * @param variants             input rows of the case (may contain duplicates)
     * @param genePhenotypeScores  per-gene phenotype score table for this case; genes absent from
     *                             the map (or mapped to {@code null}) are unscoreable;
     *                             values must be finite; may be {@code null}, treated as empty
     * @param settings             weight and missing policies
     * @return one {@link VariantScore} per distinct variant, sorted by descending combined
     *         score with a deterministic variant-key tiebreak for presentation only
     */
    public List<VariantScore> score(List<VariantInput> variants, Map<String, Double> genePhenotypeScores, Settings settings) {
        Objects.requireNonNull(variants, "variants must not be null");
        Objects.requireNonNull(settings, "settings must not be null");
        Map<String, Double> phenotypeTable = genePhenotypeScores == null ? Map.of() : new HashMap<>(genePhenotypeScores);
        for (Map.Entry<String, Double> entry : phenotypeTable.entrySet()) {
            if (entry.getValue() != null && !Double.isFinite(entry.getValue())) {
                throw new IllegalArgumentException("phenotype score for gene " + entry.getKey() + " must be finite");
            }
        }

        if (variants.isEmpty()) {
            return List.of();
        }

        // Merge duplicate rows by genomic key (insertion order kept for determinism).
        Map<VariantKey, Set<String>> genesByKey = new LinkedHashMap<>();
        Map<VariantKey, Double> pathogenicityByKey = new LinkedHashMap<>();
        for (VariantInput input : variants) {
            Objects.requireNonNull(input, "variants must not contain null entries");
            genesByKey.computeIfAbsent(input.key(), k -> new LinkedHashSet<>()).addAll(input.geneIds());
            Double gp = input.pathogenicityScore();
            Double current = pathogenicityByKey.get(input.key());
            if (gp != null && (current == null || gp > current)) {
                pathogenicityByKey.put(input.key(), gp);
            }
        }

        // Raw phenotype score per variant: max over its genes present in the table.
        Map<VariantKey, Double> rawPhenotypeByKey = new LinkedHashMap<>();
        for (Map.Entry<VariantKey, Set<String>> entry : genesByKey.entrySet()) {
            Double max = null;
            for (String geneId : entry.getValue()) {
                Double geneScore = phenotypeTable.get(geneId);
                if (geneScore == null) {
                    continue;
                }
                if (!Double.isFinite(geneScore)) {
                    throw new IllegalArgumentException(
                            "phenotype score for gene " + geneId + " must be finite, got " + geneScore);
                }
                if (max == null || geneScore > max) {
                    max = geneScore;
                }
            }
            rawPhenotypeByKey.put(entry.getKey(), max);
        }

        List<VariantKey> keys = new ArrayList<>(genesByKey.keySet());
        int n = keys.size();

        // Impute missing values per policy (median over the scoreable distinct variants).
        double pathogenicityFill = fill(settings.pathogenicityMissing(),
                pathogenicityByKey.values().stream().filter(Objects::nonNull).toList());
        double phenotypeFill = fill(settings.phenotypeMissing(),
                rawPhenotypeByKey.values().stream().filter(Objects::nonNull).toList());

        double[] imputedGp = new double[n];
        double[] imputedMs = new double[n];
        for (int i = 0; i < n; i++) {
            VariantKey key = keys.get(i);
            Double gp = pathogenicityByKey.get(key);
            Double ms = rawPhenotypeByKey.get(key);
            imputedGp[i] = gp == null ? pathogenicityFill : gp;
            imputedMs[i] = ms == null ? phenotypeFill : ms;
        }

        double[] normalizedGp = minMax(imputedGp);
        double[] normalizedMs = minMax(imputedMs);

        double[] combined = new double[n];
        double w = settings.phenotypeWeight();
        for (int i = 0; i < n; i++) {
            combined[i] = w * normalizedMs[i] + (1.0 - w) * normalizedGp[i];
            if (!Double.isFinite(combined[i])) {
                throw new IllegalArgumentException("combined score must be finite");
            }
        }

        // Exact score equality defines ties; average tie rank as in the reference harness.
        double[] rank = new double[n];
        double[] rankStrict = new double[n];
        for (int i = 0; i < n; i++) {
            int nGreater = 0;
            int nTied = 0;
            for (int j = 0; j < n; j++) {
                if (combined[j] > combined[i]) {
                    nGreater++;
                } else if (combined[j] == combined[i]) {
                    nTied++;
                }
            }
            rank[i] = nGreater + (nTied + 1) / 2.0;
            rankStrict[i] = nGreater + 1.0;
        }

        List<VariantScore> scores = new ArrayList<>(n);
        for (int i = 0; i < n; i++) {
            VariantKey key = keys.get(i);
            scores.add(new VariantScore(
                    key,
                    genesByKey.get(key),
                    rawPhenotypeByKey.get(key),
                    pathogenicityByKey.get(key),
                    rawPhenotypeByKey.get(key) == null,
                    pathogenicityByKey.get(key) == null,
                    imputedMs[i],
                    imputedGp[i],
                    normalizedMs[i],
                    normalizedGp[i],
                    combined[i],
                    rank[i],
                    rankStrict[i]));
        }

        scores.sort(Comparator.comparing(VariantScore::combinedScore, Comparator.reverseOrder())
                .thenComparing(VariantScore::key));
        return List.copyOf(scores);
    }

    private static double fill(MissingPolicy policy, List<Double> known) {
        return switch (policy) {
            case ZERO -> 0.0;
            case MEDIAN -> known.isEmpty() ? 0.0 : median(known);
        };
    }

    /**
     * Median as in Python's {@code statistics.median}: middle value for odd counts, mean of the
     * two middle values for even counts.
     */
    private static double median(List<Double> values) {
        double[] sorted = values.stream().mapToDouble(Double::doubleValue).sorted().toArray();
        int n = sorted.length;
        return n % 2 == 1 ? sorted[n / 2] : sorted[n / 2 - 1] / 2.0 + sorted[n / 2] / 2.0;
    }

    /**
     * Per-case min-max rescaling; a constant column maps to 0.
     */
    private static double[] minMax(double[] values) {
        double[] out = new double[values.length];
        if (values.length == 0) {
            return out;
        }
        double lo = values[0];
        double hi = values[0];
        for (double v : values) {
            if (v < lo) {
                lo = v;
            }
            if (v > hi) {
                hi = v;
            }
        }
        double span = hi - lo;
        if (!Double.isFinite(span)) {
            throw new IllegalArgumentException("component range must be finite after imputation");
        }
        for (int i = 0; i < values.length; i++) {
            out[i] = span == 0 ? 0.0 : (values[i] - lo) / span;
        }
        return out;
    }
}
