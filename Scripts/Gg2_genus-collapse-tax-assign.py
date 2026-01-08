#!/usr/bin/env python

import os
import pandas as pd
import biom
from biom import load_table
from Bio import SeqIO


def extract_genus_only(taxon):
    if pd.isna(taxon):
        return None

    for part in taxon.split(";"):
        part = part.strip()
        if part.startswith("g__") and part != "g__":
            return part

    return None


REGIONS = {
    "V4": {
        "biom": "../Data/16S/Tables/from_Qiita/174951_feature-table_16S_V4_rare-3769.biom",
        "fasta": "../Data/16S/Fasta/174951_V4_ASVs.fasta",
        "taxonomy": "../Analyses/16S/SEPP/174951_V4_GG2_taxonomy/taxonomy.tsv",
        "mapping_out": "../Data/16S/ASV_to_Genus_V4.tsv"
    },
    "V1-V3": {
        "biom": "../Data/16S/Tables/from_Qiita/179426_feature-table_16S_V1V3_rare-11054.biom",
        "fasta": "../Data/16S/Fasta/179426_V1V3_ASVs.fasta",
        "taxonomy": "../Analyses/16S/SEPP/179426_V1V3_GG2_taxonomy/taxonomy.tsv",
        "mapping_out": "../Data/16S/ASV_to_Genus_V1V3.tsv"
    }
}

OUT_DIR = "../Data/16S/Tables/with_taxonomy"
os.makedirs(OUT_DIR, exist_ok=True)


for region, cfg in REGIONS.items():
    print(f"\nProcessing {region}")

    table = load_table(cfg["biom"])
    obs_df = table.to_dataframe(dense=True)

    seq_to_feature = {
        str(rec.seq): rec.id
        for rec in SeqIO.parse(cfg["fasta"], "fasta")
    }

    meta_df = pd.DataFrame({"ASV_sequence": obs_df.index})
    meta_df["Feature_ID"] = meta_df["ASV_sequence"].map(seq_to_feature)

    tax_df = (
        pd.read_csv(cfg["taxonomy"], sep="\t")
        .rename(columns={"Feature ID": "Feature_ID"})
    )

    meta_df = meta_df.merge(
        tax_df[["Feature_ID", "Taxon"]],
        on="Feature_ID",
        how="left"
    )

    meta_df["Genus"] = meta_df["Taxon"].apply(extract_genus_only)
    meta_df = meta_df.dropna(subset=["Genus"]).copy()

    obs_df = obs_df.loc[meta_df["ASV_sequence"]]
    obs_df["Genus"] = meta_df.set_index("ASV_sequence")["Genus"]

    genus_df = (
        obs_df
        .groupby("Genus")
        .sum()
    )

    meta_df.to_csv(cfg["mapping_out"], sep="\t", index=False)

    genus_table = biom.Table(
        genus_df.values,
        observation_ids=genus_df.index.tolist(),
        sample_ids=genus_df.columns.tolist()
    )

    biom_out = os.path.join(
        OUT_DIR,
        os.path.basename(cfg["biom"]).replace(".biom", "_Gg2_genus_collapsed.biom")
    )

    with biom.util.biom_open(biom_out, "w") as f:
        genus_table.to_hdf5(
            f,
            "Greengenes2 genus-collapsed table"
        )

    print(f"  Genera retained: {genus_table.shape[0]}")
    print(f"  Saved BIOM: {biom_out}")

print("\nDone.")
