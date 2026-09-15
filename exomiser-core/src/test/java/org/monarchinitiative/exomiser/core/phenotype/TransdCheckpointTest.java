package org.monarchinitiative.exomiser.core.phenotype;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;

class TransdCheckpointTest {

    @TempDir Path directory;

    @Test
    void scoresProjectionAndSquaredDistance() throws IOException {
        fixture("0.5,0.0", "0.0,2.0", "0.0,0.0", "0.0,0.0",
                "0.0,0.0", "0.0,1.0");
        // Head: e=(.5,0), dot(e_p,e)=0, projected=(.5,0).
        // Relation is zero, tail zero, so raw -L2^2 = -.25.
        assertEquals(-0.25, TransdCheckpoint.load(directory).scoreGeneCase("123", "PAVS:A1"), 1e-12);
    }

    @Test
    void clampsProjectedEntityToUnitNorm() throws IOException {
        fixture("1.0,0.0", "2.0,0.0", "0.0,0.0", "0.0,0.0",
                "0.0,0.0", "1.0,0.0");
        // Head projection=(3,0) before clamp, then (1,0).
        assertEquals(-1.0, TransdCheckpoint.load(directory).scoreGeneCase("123", "PAVS:A1"), 1e-12);
    }

    @Test
    void usesRelationOffsetAndTailProjection() throws IOException {
        fixture("0.5,0.0", "0.0,0.0", "0.0,0.5", "0.0,0.5",
                "0.25,0.0", "0.0,1.0");
        // Head=(.5,0), tail=(0,.75), relation=(.25,0).
        // Residual=(.75,-.75), score=-(.75^2+.75^2).
        assertEquals(-1.125, TransdCheckpoint.load(directory).scoreGeneCase("123", "PAVS:A1"), 1e-12);
    }

    @Test
    void missingCaseFailsClearly() throws IOException {
        fixture("1.0,0.0", "0.0,0.0", "0.0,0.0", "0.0,0.0",
                "0.0,0.0", "0.0,0.0");
        assertThrows(IllegalArgumentException.class,
                () -> TransdCheckpoint.load(directory).scoreGeneCase("123", "PAVS:UNKNOWN"));
    }

    @Test
    void reportsVocabularyPresenceSeparately() throws IOException {
        fixture("1.0,0.0", "0.0,0.0", "0.0,0.0", "0.0,0.0",
                "0.0,0.0", "0.0,0.0");
        TransdCheckpoint checkpoint = TransdCheckpoint.load(directory);
        assertEquals(true, checkpoint.hasGene("123"));
        assertEquals(false, checkpoint.hasGene("999"));
        assertEquals(true, checkpoint.hasCase("PAVS:A1"));
        assertEquals(false, checkpoint.hasCase("PAVS:UNKNOWN"));
    }

    private void fixture(String gene, String geneProjection, String patient,
                         String patientProjection, String relation, String relationProjection)
            throws IOException {
        Files.writeString(directory.resolve("entities.tsv"),
                "http://mowl.borg/123\t" + gene + "\t" + geneProjection + "\n" +
                "http://mowl.borg/PAVS_A1\t" + patient + "\t" + patientProjection + "\n");
        Files.writeString(directory.resolve("relations.tsv"),
                "http://mowl.borg/associated_with\t" + relation + "\t" + relationProjection + "\n");
        Files.writeString(directory.resolve("manifest.properties"),
                "format=pykeen-transd-v1\npykeen_version=1.11.1\n" +
                "model_state_sha256=" + "0".repeat(64) + "\n" +
                "entity_dimension=2\nrelation_dimension=2\nentity_count=2\nrelation_count=1\n" +
                "entities_sha256=" + sha256(directory.resolve("entities.tsv")) + "\n" +
                "relations_sha256=" + sha256(directory.resolve("relations.tsv")) + "\n" +
                "norm_p=2\npower_norm=true\nprojection_maxnorm=1\n" +
                "predict_with_sigmoid=false\n" +
                "association_relation=http://mowl.borg/associated_with\n");
    }

    private static String sha256(Path file) throws IOException {
        try {
            byte[] bytes = MessageDigest.getInstance("SHA-256").digest(Files.readAllBytes(file));
            StringBuilder hex = new StringBuilder(64);
            for (byte b : bytes) {
                hex.append(Character.forDigit((b >>> 4) & 0xf, 16));
                hex.append(Character.forDigit(b & 0xf, 16));
            }
            return hex.toString();
        } catch (NoSuchAlgorithmException e) {
            throw new IOException(e);
        }
    }
}
