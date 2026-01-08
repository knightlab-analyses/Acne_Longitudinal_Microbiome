#!/usr/bin/env python

import os
import pandas as pd
import biom
from biom import load_table
from Bio import SeqIO


def extract_highest_taxon(taxon):
    if pd.isna(taxon):
        return None

    ranks = ["g__", "f__", "o__", "c__", "p__", "k__"]

    parts = [p.strip() for p in taxon.split(";")]

    for rank in ranks:
        for part in parts:
            if part.startswith(rank) and part != rank:
                return part

    return None


REGIONS = {
    "V4": {
        "biom": "../Data/16S/Tables/from_Qiita/174951_feature-table_16S_V4_rare-3769.biom",
        "fasta": "../Data/16S/Fasta/174951_V4_ASVs.fasta",
        "taxonomy": "../Analyses/16S/SEPP/174951_V4_GG2_taxonomy/taxonomy.tsv",
        "mapping_out": "../Analyses/16S/SEPP/174951_V4_GG2_taxonomy/ASV_to_Taxon_V4.tsv"
    },
    "V1-V3": {
        "biom": "../Data/16S/Tables/from_Qiita/179426_feature-table_16S_V1V3_rare-11054.biom",
        "fasta": "../Data/16S/Fasta/179426_V1V3_ASVs.fasta",
        "taxonomy": "../Analyses/16S/SEPP/179426_V1V3_GG2_taxonomy/taxonomy.tsv",
        "mapping_out": "../Analyses/16S/SEPP/179426_V1V3_GG2_taxonomy/ASV_to_Taxon_V1V3.tsv"
    }
}

OUT_DIR = "../Data/16S/Tables/with_taxonomy"
os.makedirs(OUT_DIR, exist_ok=True)


for region, cfg in REGIONS.items():
    print(f"\nProcessing {region}")

    table = load_table(cfg["biom"])
    obs_df = table.to_dataframe(dense=True)

    seq_to_feature = {
        str(rec.seq).upper(): rec.id
        for rec in SeqIO.parse(cfg["fasta"], "fasta")
    }

    meta_df = pd.DataFrame({"ASV_sequence": obs_df.index.str.upper()})
    meta_df["Feature_ID"] = meta_df["ASV_sequence"].map(seq_to_feature)

    missing_fasta = meta_df["Feature_ID"].isna().sum()
    if missing_fasta:
        print(f"  Warning: {missing_fasta} ASVs not found in FASTA")

    tax_df = (
        pd.read_csv(cfg["taxonomy"], sep="\t")
        .rename(columns={"Feature ID": "Feature_ID"})
    )

    meta_df = meta_df.merge(
        tax_df[["Feature_ID", "Taxon"]],
        on="Feature_ID",
        how="left"
    )

    meta_df["Taxon_label"] = meta_df["Taxon"].apply(extract_highest_taxon)

    before = meta_df.shape[0]
    meta_df = meta_df.dropna(subset=["Taxon_label"])
    after = meta_df.shape[0]
    print(f"  Retained {after}/{before} ASVs with any taxonomic assignment")

    obs_df = obs_df.loc[meta_df["ASV_sequence"]]

    obs_df["Taxon_label"] = meta_df.set_index("ASV_sequence")["Taxon_label"]
    obs_df["total_abundance"] = obs_df.drop(columns="Taxon_label").sum(axis=1)

    obs_df = obs_df.sort_values(
        by=["Taxon_label", "total_abundance"],
        ascending=[True, False]
    )

    obs_df["ASV_rank"] = obs_df.groupby("Taxon_label").cumcount().add(1)

    obs_df["New_Feature_ID"] = (
        obs_df["Taxon_label"] + "_ASV-" + obs_df["ASV_rank"].astype(str)
    )

    mapping_df = obs_df.reset_index()[[
        "index", "Taxon_label", "ASV_rank", "total_abundance"
    ]].rename(columns={"index": "ASV_sequence"})

    mapping_df.to_csv(cfg["mapping_out"], sep="\t", index=False)

    final_df = obs_df.drop(
        columns=["Taxon_label", "ASV_rank", "total_abundance", "New_Feature_ID"]
    )
    final_df.index = obs_df["New_Feature_ID"]

    renamed_table = biom.Table(
        final_df.values,
        observation_ids=final_df.index.tolist(),
        sample_ids=final_df.columns.tolist()
    )

    biom_out = os.path.join(
        OUT_DIR,
        os.path.basename(cfg["biom"]).replace(".biom", "_Gg2_highest-taxon-ASV.biom")
    )

    with biom.util.biom_open(biom_out, "w") as f:
        renamed_table.to_hdf5(
            f,
            "Greengenes2 taxon-ranked ASVs (highest available rank)"
        )

    print(f"  Final features: {renamed_table.shape[0]}")
    print(f"  Saved BIOM: {biom_out}")

print("\nDone.")
