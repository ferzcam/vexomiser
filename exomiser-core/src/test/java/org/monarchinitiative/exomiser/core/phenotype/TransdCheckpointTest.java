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
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.junit.jupiter.api.Assumptions.assumeTrue;

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

    @Test
    void rejectsInverseTripleManifest() throws IOException {
        fixture("1.0,0.0", "0.0,0.0", "0.0,0.0", "0.0,0.0",
                "0.0,0.0", "0.0,0.0");
        Path manifest = directory.resolve("manifest.properties");
        Files.writeString(manifest, Files.readString(manifest).replace(
                "inverse_triples=false", "inverse_triples=true"));
        assertThrows(IOException.class, () -> TransdCheckpoint.load(directory));
    }

    @Test
    void rejectsModifiedWeightsByChecksum() throws IOException {
        fixture("1.0,0.0", "0.0,0.0", "0.0,0.0", "0.0,0.0",
                "0.0,0.0", "0.0,0.0");
        Path entities = directory.resolve("entities.tsv");
        Files.writeString(entities, Files.readString(entities).replace("1.0,0.0", "0.5,0.0"));
        assertThrows(IOException.class, () -> TransdCheckpoint.load(directory));
    }

    @Test
    void matchesSavedPykeenScoresWhenWorkstationFixtureIsSupplied() throws IOException {
        String bundle = System.getProperty("transd.parity.bundle");
        String expectedFile = System.getProperty("transd.parity.expected");
        String actualFile = System.getProperty("transd.parity.actual");
        assumeTrue(bundle != null || expectedFile != null || actualFile != null,
                "workstation parity fixture paths were not supplied");
        assertTrue(bundle != null && expectedFile != null && actualFile != null,
                "parity requires bundle, expected, and actual output paths");
        assertTrue(Files.exists(Path.of(bundle)) && Files.exists(Path.of(expectedFile)),
                "parity bundle and expected scores must exist");
        TransdCheckpoint checkpoint = TransdCheckpoint.load(Path.of(bundle));
        var lines = Files.readAllLines(Path.of(expectedFile));
        if (lines.isEmpty() || !"gene_iri\tcase_iri\traw_pykeen_score".equals(lines.getFirst())) {
            throw new IOException("invalid PyKEEN parity expected header");
        }
        StringBuilder result = new StringBuilder("gene_iri\tcase_iri\traw_pykeen_score\traw_java_score\tabs_error\n");
        double maxAbsoluteError = 0;
        for (String line : lines.subList(1, lines.size())) {
            String[] fields = line.split("\\t", -1);
            if (fields.length != 3) {
                throw new IOException("invalid PyKEEN parity expected row: " + line);
            }
            String geneId = fields[0].substring("http://mowl.borg/".length());
            String caseId = fields[1].substring("http://mowl.borg/".length()).replace('_', ':');
            double expected = Double.parseDouble(fields[2]);
            double actual = checkpoint.scoreGeneCase(geneId, caseId);
            double error = Math.abs(expected - actual);
            maxAbsoluteError = Math.max(maxAbsoluteError, error);
            result.append(fields[0]).append('\t').append(fields[1]).append('\t')
                    .append(expected).append('\t').append(actual).append('\t').append(error).append('\n');
        }
        Path actualPath = Path.of(actualFile);
        Files.createDirectories(actualPath.toAbsolutePath().getParent());
        Files.writeString(actualPath, result.toString());
        System.out.println("TransD PyKEEN/Java parity rows=" + (lines.size() - 1)
                + " max absolute error=" + maxAbsoluteError + " actual=" + actualPath);
        assertTrue(lines.size() > 1, "PyKEEN parity expected file has no prediction rows");
        assertTrue(maxAbsoluteError <= 1e-5,
                "PyKEEN/Java TransD parity exceeded 1e-5; inspect saved actual scores");
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
                "norm_p=2\npower_norm=true\nprojection_maxnorm=1\ninverse_triples=false\n" +
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
