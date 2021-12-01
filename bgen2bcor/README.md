# Generate BCOR files from a bgen file

This WDL generates per-chromosome BCOR files from a bgen file. LDstore v1.1 is used.

## WDL steps

1. Bgen file is split to chromosomes with bgenix and indexed.
2. Per-chromosome bgens are used to generate per-chromosome BCOR files.

The example config json was used to generate BCORs from the SISuv4 imputation panel.

## Usage

In the config, define your bgen file to use (`bgen2bcor.bgen`).
If needed, adjust the LDstore parameters `variant_window_size`, `ld_thold`, and `low` to your liking.
