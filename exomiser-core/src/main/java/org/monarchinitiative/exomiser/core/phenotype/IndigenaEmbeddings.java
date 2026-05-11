package org.monarchinitiative.exomiser.core.phenotype;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.HashMap;
import java.util.Map;

/**
 * Loads pre-computed INDIGENA entity embeddings exported by eval/export_embeddings.py.
 *
 * TSV format: {@code <entity_IRI>\t<v1>,<v2>,...,<vD>}
 *
 * Phenotype IDs from Exomiser (e.g. {@code HP:0001234}, {@code MP:0000633}) are
 * converted to OBO IRIs ({@code http://purl.obolibrary.org/obo/HP_0001234}) before lookup.
 */
public class IndigenaEmbeddings {

    private final Map<String, float[]> store;
    private final int dim;

    private IndigenaEmbeddings(Map<String, float[]> store, int dim) {
        this.store = store;
        this.dim = dim;
    }

    public static IndigenaEmbeddings load(Path tsvPath) throws IOException {
        Map<String, float[]> map = new HashMap<>();
        int dim = -1;
        try (BufferedReader br = Files.newBufferedReader(tsvPath)) {
            String line;
            while ((line = br.readLine()) != null) {
                int tab = line.indexOf('\t');
                if (tab < 0) continue;
                String iri = line.substring(0, tab);
                String[] parts = line.substring(tab + 1).split(",");
                float[] vec = new float[parts.length];
                for (int i = 0; i < parts.length; i++) {
                    vec[i] = Float.parseFloat(parts[i]);
                }
                if (dim < 0) dim = vec.length;
                map.put(iri, vec);
            }
        }
        return new IndigenaEmbeddings(map, Math.max(dim, 0));
    }

    /**
     * Convert a phenotype ID in Exomiser format to an OBO IRI.
     * e.g. {@code HP:0001234} → {@code http://purl.obolibrary.org/obo/HP_0001234}
     */
    public static String toOboIri(String phenotypeId) {
        return "http://purl.obolibrary.org/obo/" + phenotypeId.replace(':', '_');
    }

    /**
     * Look up the embedding for a phenotype given its Exomiser ID (e.g. {@code HP:0001234}).
     * Returns {@code null} if the term is not in the embedding vocabulary.
     */
    public float[] getByPhenotypeId(String phenotypeId) {
        return store.get(toOboIri(phenotypeId));
    }

    /** Look up the embedding by full OBO IRI. Returns {@code null} if absent. */
    public float[] getByIri(String iri) {
        return store.get(iri);
    }

    public int dim() {
        return dim;
    }

    public int size() {
        return store.size();
    }
}
