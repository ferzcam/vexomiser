package org.monarchinitiative.exomiser.core.phenotype;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedHashSet;
import java.util.Map;
import java.util.Set;

/**
 * Gene and case embeddings for the EmbedPVP phenotype score.
 *
 * <p>EmbedPVP's phenotype term is an <em>entity</em>-level similarity — the case embedding
 * against the gene embedding — not a comparison of phenotype term sets. The case
 * embedding is produced by an offline inductive fine-tune and so must be exported and
 * read back here; see {@code docs/embedpvp_scorer.md} for the contract this class
 * implements.
 *
 * <p>Gene file: {@code http://mowl.borg/<entrezGeneId>\t<v1>,<v2>,...}
 *
 * <p>Case file: {@code <caseId>\t<v1>,<v2>,...[\t<HP:...;HP:...>]} where {@code caseId} is
 * the case id as written in the phenopacket (e.g. {@code PAVS:A0000003}) and the optional
 * third column is the sorted HPO set the embedding was computed from, used only as a
 * wiring check.
 */
public class EmbedpvpEmbeddings {

    /** IRI prefix mOWL gives graph entities; genes are keyed by Entrez id under it. */
    public static final String GENE_IRI_PREFIX = "http://mowl.borg/";

    private final Map<Integer, float[]> geneVectors;
    private final Map<String, float[]> caseVectors;
    private final Map<String, Set<String>> casePhenotypes;
    private final int dim;

    private EmbedpvpEmbeddings(Map<Integer, float[]> geneVectors,
                               Map<String, float[]> caseVectors,
                               Map<String, Set<String>> casePhenotypes,
                               int dim) {
        this.geneVectors = geneVectors;
        this.caseVectors = caseVectors;
        this.casePhenotypes = casePhenotypes;
        this.dim = dim;
    }

    /**
     * @param geneEmbeddingsTsv gene embeddings, keyed by {@code http://mowl.borg/<entrez>}
     * @param caseEmbeddingsTsv case embeddings, keyed by case id; may be null when the
     *                          gene file also carries the case rows (transductive models)
     */
    public static EmbedpvpEmbeddings load(Path geneEmbeddingsTsv, Path caseEmbeddingsTsv)
            throws IOException {
        Map<Integer, float[]> genes = new HashMap<>();
        Map<String, float[]> cases = new HashMap<>();
        Map<String, Set<String>> casePhenos = new HashMap<>();
        int[] dim = {-1};

        readGeneFile(geneEmbeddingsTsv, genes, cases, dim);
        if (caseEmbeddingsTsv != null) {
            readCaseFile(caseEmbeddingsTsv, cases, casePhenos, dim);
        }
        return new EmbedpvpEmbeddings(genes, cases, casePhenos, Math.max(dim[0], 0));
    }

    private static void readGeneFile(Path path, Map<Integer, float[]> genes,
                                     Map<String, float[]> cases, int[] dim)
            throws IOException {
        try (BufferedReader br = Files.newBufferedReader(path)) {
            String line;
            while ((line = br.readLine()) != null) {
                int tab = line.indexOf('\t');
                if (tab < 0) {
                    continue;
                }
                String key = line.substring(0, tab);
                if (!key.startsWith(GENE_IRI_PREFIX)) {
                    continue;
                }
                String local = key.substring(GENE_IRI_PREFIX.length());
                float[] vec = parseVector(line, tab + 1, dim);
                Integer entrez = parseEntrez(local);
                if (entrez != null) {
                    genes.put(entrez, vec);
                } else {
                    // Transductive exports carry case entities under the same prefix, with
                    // the colon mangled to an underscore (PAVS_A0000003). Restore it so the
                    // same file can serve as the case source.
                    cases.putIfAbsent(local.replace('_', ':'), vec);
                }
            }
        }
    }

    private static void readCaseFile(Path path, Map<String, float[]> cases,
                                     Map<String, Set<String>> casePhenos, int[] dim)
            throws IOException {
        try (BufferedReader br = Files.newBufferedReader(path)) {
            String line;
            while ((line = br.readLine()) != null) {
                String[] fields = line.split("\t");
                if (fields.length < 2) {
                    continue;
                }
                float[] vec = parseVector(fields[1], 0, dim);
                cases.put(fields[0], vec);           // explicit file wins over the gene file
                if (fields.length >= 3 && !fields[2].isBlank()) {
                    Set<String> terms = new LinkedHashSet<>();
                    for (String t : fields[2].split(";")) {
                        if (!t.isBlank()) {
                            terms.add(t.trim());
                        }
                    }
                    casePhenos.put(fields[0], terms);
                }
            }
        }
    }

    private static float[] parseVector(String source, int from, int[] dim) {
        String[] parts = source.substring(from).trim().split(",");
        float[] vec = new float[parts.length];
        for (int i = 0; i < parts.length; i++) {
            vec[i] = Float.parseFloat(parts[i]);
        }
        if (dim[0] < 0) {
            dim[0] = vec.length;
        }
        return vec;
    }

    private static Integer parseEntrez(String local) {
        for (int i = 0; i < local.length(); i++) {
            if (local.charAt(i) < '0' || local.charAt(i) > '9') {
                return null;
            }
        }
        return local.isEmpty() ? null : Integer.valueOf(local);
    }

    /** Gene embedding for an Entrez gene id, or null if the gene is not in the vocabulary. */
    public float[] getGeneVector(int entrezGeneId) {
        return geneVectors.get(entrezGeneId);
    }

    /** Case embedding for a case id (e.g. {@code PAVS:A0000003}), or null if absent. */
    public float[] getCaseVector(String caseId) {
        return caseId == null ? null : caseVectors.get(caseId);
    }

    /**
     * The HPO term set recorded for a case at export time, empty when the export did not
     * carry one. Used only to detect a case embedding wired to the wrong patient.
     */
    public Set<String> getCasePhenotypes(String caseId) {
        Set<String> terms = casePhenotypes.get(caseId);
        return terms == null ? Collections.emptySet() : terms;
    }

    public int dim() {
        return dim;
    }

    public int geneCount() {
        return geneVectors.size();
    }

    public int caseCount() {
        return caseVectors.size();
    }
}
