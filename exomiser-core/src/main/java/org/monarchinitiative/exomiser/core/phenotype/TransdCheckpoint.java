package org.monarchinitiative.exomiser.core.phenotype;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.HashMap;
import java.util.Map;
import java.util.Properties;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;

/** Exact raw PyKEEN TransD triple inference from a selected checkpoint bundle. */
public final class TransdCheckpoint {

    private record Pair(double[] embedding, double[] projection) {}

    private final Map<String, Pair> entities;
    private final Pair association;
    private final int relationDimension;

    private TransdCheckpoint(Map<String, Pair> entities, Pair association) {
        this.entities = Map.copyOf(entities);
        this.association = association;
        this.relationDimension = association.embedding().length;
    }

    /** Load the files written by code/gda/export_transd_bundle.py. */
    public static TransdCheckpoint load(Path directory) throws IOException {
        Properties manifest = new Properties();
        try (var reader = Files.newBufferedReader(directory.resolve("manifest.properties"))) {
            manifest.load(reader);
        }
        require(manifest, "format", "pykeen-transd-v1");
        require(manifest, "norm_p", "2");
        require(manifest, "power_norm", "true");
        require(manifest, "projection_maxnorm", "1");
        require(manifest, "predict_with_sigmoid", "false");
        int entityDimension = positiveInt(manifest, "entity_dimension");
        int relationDimension = positiveInt(manifest, "relation_dimension");
        int entityCount = positiveInt(manifest, "entity_count");
        int relationCount = positiveInt(manifest, "relation_count");
        if (manifest.getProperty("pykeen_version") == null) {
            throw new IOException("TransD bundle lacks PyKEEN version provenance");
        }
        if (manifest.getProperty("model_state_sha256") == null ||
                !manifest.getProperty("model_state_sha256").matches("[a-f0-9]{64}")) {
            throw new IOException("TransD bundle lacks model state provenance");
        }
        String relationIri = manifest.getProperty("association_relation");
        if (relationIri == null || relationIri.isBlank()) {
            throw new IOException("TransD bundle lacks association_relation");
        }
        Path entityFile = directory.resolve("entities.tsv");
        Path relationFile = directory.resolve("relations.tsv");
        verifyChecksum(entityFile, manifest.getProperty("entities_sha256"));
        verifyChecksum(relationFile, manifest.getProperty("relations_sha256"));
        Map<String, Pair> entities = readRows(entityFile, entityDimension);
        Map<String, Pair> relations = readRows(relationFile, relationDimension);
        if (entities.size() != entityCount || relations.size() != relationCount) {
            throw new IOException("TransD bundle row counts disagree with manifest");
        }
        Pair association = relations.get(relationIri);
        if (association == null) {
            throw new IOException("TransD association relation is absent: " + relationIri);
        }
        return new TransdCheckpoint(entities, association);
    }

    /** Raw score for (gene, associated_with, case); higher is better. */
    public double scoreGeneCase(String entrezGeneId, String caseId) {
        if (entrezGeneId == null || entrezGeneId.isBlank() || caseId == null || caseId.isBlank()) {
            throw new IllegalArgumentException("TransD gene and case IDs are required");
        }
        Pair gene = entity("http://mowl.borg/" + entrezGeneId);
        Pair patient = entity("http://mowl.borg/" + caseId.replace(':', '_'));
        double[] head = project(gene);
        double[] tail = project(patient);
        double squaredDistance = 0;
        for (int i = 0; i < relationDimension; i++) {
            double delta = head[i] + association.embedding()[i] - tail[i];
            squaredDistance += delta * delta;
        }
        return -squaredDistance;
    }

    public boolean hasGene(String entrezGeneId) {
        return entrezGeneId != null && entities.containsKey("http://mowl.borg/" + entrezGeneId);
    }

    public boolean hasCase(String caseId) {
        return caseId != null && entities.containsKey("http://mowl.borg/" + caseId.replace(':', '_'));
    }

    private Pair entity(String iri) {
        Pair value = entities.get(iri);
        if (value == null) {
            throw new IllegalArgumentException("TransD checkpoint has no entity: " + iri);
        }
        return value;
    }

    // PyKEEN utils.project_entity: r_p * dot(e_p,e) + identity-padded e,
    // then clamp_norm(p=2, maxnorm=1). No normalization of stored weights here.
    private double[] project(Pair entity) {
        double dot = 0;
        for (int i = 0; i < entity.embedding().length; i++) {
            dot += entity.projection()[i] * entity.embedding()[i];
        }
        double[] result = new double[relationDimension];
        double squaredNorm = 0;
        for (int i = 0; i < relationDimension; i++) {
            result[i] = association.projection()[i] * dot;
            if (i < entity.embedding().length) {
                result[i] += entity.embedding()[i];
            }
            squaredNorm += result[i] * result[i];
        }
        if (squaredNorm > 1) {
            double scale = 1 / Math.sqrt(squaredNorm);
            for (int i = 0; i < result.length; i++) {
                result[i] *= scale;
            }
        }
        return result;
    }

    private static Map<String, Pair> readRows(Path file, int dimension) throws IOException {
        Map<String, Pair> rows = new HashMap<>();
        int lineNumber = 0;
        try (var reader = Files.newBufferedReader(file)) {
            String line;
            while ((line = reader.readLine()) != null) {
                lineNumber++;
                String[] fields = line.split("\\t", -1);
                if (fields.length != 3 || fields[0].isBlank()) {
                    throw new IOException(file + ":" + lineNumber + " invalid TransD row");
                }
                Pair value = new Pair(vector(fields[1], dimension, file, lineNumber),
                                      vector(fields[2], dimension, file, lineNumber));
                if (rows.putIfAbsent(fields[0], value) != null) {
                    throw new IOException(file + ":" + lineNumber + " duplicate IRI: " + fields[0]);
                }
            }
        }
        if (rows.isEmpty()) {
            throw new IOException("TransD bundle has no rows: " + file);
        }
        return rows;
    }

    private static double[] vector(String text, int dimension, Path file, int lineNumber) throws IOException {
        String[] parts = text.split(",", -1);
        if (parts.length != dimension) {
            throw new IOException(file + ":" + lineNumber + " expected " + dimension + " components");
        }
        double[] values = new double[dimension];
        for (int i = 0; i < dimension; i++) {
            try {
                values[i] = Double.parseDouble(parts[i]);
            } catch (NumberFormatException e) {
                throw new IOException(file + ":" + lineNumber + " invalid component", e);
            }
            if (!Double.isFinite(values[i])) {
                throw new IOException(file + ":" + lineNumber + " non-finite component");
            }
        }
        return values;
    }

    private static void require(Properties properties, String key, String expected) throws IOException {
        if (!expected.equals(properties.getProperty(key))) {
            throw new IOException("TransD bundle requires " + key + "=" + expected);
        }
    }

    private static int positiveInt(Properties properties, String key) throws IOException {
        try {
            int value = Integer.parseInt(properties.getProperty(key));
            if (value > 0) return value;
        } catch (NumberFormatException ignored) {
            // Report one consistent error for missing and malformed dimensions.
        }
        throw new IOException("TransD bundle requires positive " + key);
    }

    private static void verifyChecksum(Path file, String expected) throws IOException {
        if (expected == null || !expected.matches("[a-f0-9]{64}")) {
            throw new IOException("TransD bundle lacks valid checksum for " + file.getFileName());
        }
        try {
            MessageDigest digest = MessageDigest.getInstance("SHA-256");
            try (var input = Files.newInputStream(file)) {
                byte[] buffer = new byte[8192];
                int length;
                while ((length = input.read(buffer)) != -1) {
                    digest.update(buffer, 0, length);
                }
            }
            StringBuilder hex = new StringBuilder(64);
            for (byte b : digest.digest()) {
                hex.append(Character.forDigit((b >>> 4) & 0xf, 16));
                hex.append(Character.forDigit(b & 0xf, 16));
            }
            if (!expected.equals(hex.toString())) {
                throw new IOException("TransD bundle checksum mismatch for " + file.getFileName());
            }
        } catch (NoSuchAlgorithmException e) {
            throw new IOException("SHA-256 unavailable", e);
        }
    }
}
