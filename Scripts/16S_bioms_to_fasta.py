#!/usr/bin/env python

import os
import hashlib
from biom import load_table

# Input BIOMs and output FASTAs
biom_to_fasta = {
    "../Data/16S/Tables/from_Qiita/179426_feature-table_16S_V1V3_rare-11054.biom":
        "../Data/16S/Fasta/179426_V1V3_ASVs.fasta",
    "../Data/16S/Tables/from_Qiita/174951_feature-table_16S_V4_rare-3769.biom":
        "../Data/16S/Fasta/174951_V4_ASVs.fasta",
}

OUT_DIR = "../Data/16S/Fasta"
os.makedirs(OUT_DIR, exist_ok=True)

# Convert each BIOM to FASTA
for biom_path, fasta_out in biom_to_fasta.items():
    print(f"Processing {biom_path}")

    table = load_table(biom_path)
    asv_seqs = table.ids(axis="observation")

    with open(fasta_out, "w") as f:
        for seq in asv_seqs:
            seq = seq.upper()
            asv_hash = hashlib.md5(seq.encode()).hexdigest()[:12]
            f.write(f">ASV_{asv_hash}\n")
            f.write(f"{seq}\n")

    print(f"  Wrote {len(asv_seqs)} ASVs → {fasta_out}")