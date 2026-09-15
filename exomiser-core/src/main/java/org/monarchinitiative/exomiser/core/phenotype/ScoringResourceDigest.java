package org.monarchinitiative.exomiser.core.phenotype;

import java.io.IOException;
import java.io.InputStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.util.HexFormat;

/** SHA-256 identities of selected offline phenotype resources. */
public final class ScoringResourceDigest {
    private ScoringResourceDigest() {}

    public static String sha256(Path file) {
        if (file == null) return "";
        try {
            MessageDigest digest = MessageDigest.getInstance("SHA-256");
            try (InputStream input = Files.newInputStream(file)) {
                byte[] buffer = new byte[65536];
                int n;
                while ((n = input.read(buffer)) != -1) digest.update(buffer, 0, n);
            }
            return HexFormat.of().formatHex(digest.digest());
        } catch (IOException | NoSuchAlgorithmException e) {
            throw new IllegalArgumentException("Unable to hash scoring resource: " + file, e);
        }
    }

    public static String transdBundle(Path directory) {
        if (directory == null) return "";
        // The manifest is part of the bundle contract; retain file digests as well for verification.
        return "manifest=" + sha256(directory.resolve("manifest.properties")) + ";entities=" +
                sha256(directory.resolve("entities.tsv")) + ";relations=" + sha256(directory.resolve("relations.tsv"));
    }
}
