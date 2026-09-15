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
        int lineNumber = 0;
        try (BufferedReader br = Files.newBufferedReader(tsvPath)) {
            String line;
            while ((line = br.readLine()) != null) {
                lineNumber++;
                int tab = line.indexOf('\t');
                if (tab < 0 || tab == 0) throw new IOException(tsvPath + ":" + lineNumber + " invalid embedding row");
                String iri = line.substring(0, tab);
                String[] parts = line.substring(tab + 1).split(",");
                if (parts.length == 0 || (dim >= 0 && parts.length != dim))
                    throw new IOException(tsvPath + ":" + lineNumber + " inconsistent embedding dimension");
                float[] vec = new float[parts.length];
                for (int i = 0; i < parts.length; i++) {
                    try { vec[i] = Float.parseFloat(parts[i]); }
                    catch (NumberFormatException e) { throw new IOException(tsvPath + ":" + lineNumber + " invalid embedding", e); }
                    if (!Float.isFinite(vec[i])) throw new IOException(tsvPath + ":" + lineNumber + " non-finite embedding");
                }
                if (dim < 0) dim = vec.length;
                if (map.putIfAbsent(iri, vec) != null)
                    throw new IOException(tsvPath + ":" + lineNumber + " duplicate embedding IRI: " + iri);
            }
        }
        if (map.isEmpty()) throw new IOException("No INDIGENA embeddings in " + tsvPath);
        return new IndigenaEmbeddings(map, dim);
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
