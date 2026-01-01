#!/usr/bin/env python

import os
import pandas as pd
import biom
from biom import load_table
from biom.util import biom_open
import warnings
import logging

######################################################
# MAP ASVs TO GENUS-ASV FEATURES USING SILVA TAXONOMY
######################################################

warnings.simplefilter(action="ignore", category=FutureWarning)

# Logging
os.makedirs("../Logs", exist_ok=True)
logging.basicConfig(
    filename="../Logs/SILVA_tax-assign.log",
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s"
)

# Load SILVA taxonomy (TSV, not QZA)
logging.info("Loading SILVA taxonomy TSV...")

silva_taxonomy_path = "../Taxonomy/174116_taxonomy.tsv"
silva_taxonomy = pd.read_csv(silva_taxonomy_path, sep="\t")

# Handle common SILVA formats
if "Feature ID" in silva_taxonomy.columns:
    silva_taxonomy = silva_taxonomy.set_index("Feature ID")

if "Taxon" not in silva_taxonomy.columns:
    raise ValueError("SILVA taxonomy file must contain a 'Taxon' column.")

logging.info(
    f"SILVA taxonomy loaded with {silva_taxonomy.shape[0]} features."
)

# Load metadata (kept for consistency; not directly used here)
metadata_path = "../Metadata/metadata_final_22102024.tsv"
metadata = pd.read_csv(metadata_path, sep="\t").set_index("SampleID")

# Genus-ASV collapsing
def add_unique_tax_labels(biom_path: str, level: str = "Genus") -> pd.DataFrame:
    """
    Convert ASVs into Genus-ASV features ranked by total abundance.

    Output:
        samples × Genus-ASVs (counts)
    """

    biom_table = load_table(biom_path)
    df = biom_table.to_dataframe().T  # samples × ASVs
    logging.info(f"BIOM shape: {df.shape}")

    # Merge SILVA taxonomy
    df_tax = (
        df.T
        .merge(
            silva_taxonomy["Taxon"],
            how="left",
            left_index=True,
            right_index=True
        )
    )

    taxonomy_levels = [
        "Kingdom", "Phylum", "Class",
        "Order", "Family", "Genus", "Species"
    ]

    df_tax["Taxon"] = df_tax["Taxon"].fillna("Unknown")

    # Split SILVA taxonomy string
    tax_split = df_tax["Taxon"].str.split(";", expand=True)
    tax_split = tax_split.reindex(
        columns=range(len(taxonomy_levels)),
        fill_value=None
    )
    tax_split.columns = taxonomy_levels

    df_tax[taxonomy_levels] = tax_split

    # Drop features without genus annotation
    df_tax = df_tax[df_tax[level].notnull()]

    # Total abundance across all samples
    df_tax["TotalAbundance"] = df_tax.drop(
        columns=taxonomy_levels + ["Taxon"]
    ).sum(axis=1)

    # Rank ASVs within each genus
    index_map = {}
    for taxon, group in df_tax.groupby(level):
        group = group.sort_values("TotalAbundance", ascending=False)

        clean_taxon = (
            taxon.strip()
            .replace(":", "_")
        )

        for i, idx in enumerate(group.index, start=1):
            index_map[idx] = f"{clean_taxon} ASV-{i}"

    # Rename ASVs → Genus-ASV IDs
    df_tax = df_tax.rename(index=index_map)

    # Final counts table
    df_out = df_tax.drop(
        columns=taxonomy_levels + ["Taxon", "TotalAbundance"]
    )

    return df_out.T  # samples × Genus-ASVs


# Save BIOM
def save_biom_table(df_counts: pd.DataFrame, output_path: str):

    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    logging.info(f"BIOM shape: {df_counts.shape}")

    biom_out = biom.table.Table(
        df_counts.T.values,
        observation_ids=df_counts.columns,
        sample_ids=df_counts.index
    )

    with biom_open(output_path, "w") as f:
        biom_out.to_hdf5(
            f,
            generated_by="SILVA Genus-ASV mapping"
        )

    logging.info(f"Saved BIOM: {output_path}")


# Main
if __name__ == "__main__":

    try:
        biom_jobs = {
            "V4": {
                "input": (
                    "../Data/16S/Tables/"
                    "174951_feature-table_16S_V4_rare-3769_sampleIDfixed.biom"
                ),
                "output": (
                    "../Data/16S/Tables/"
                    "174951_feature-table_16S_V4_rare-3769_sampleIDfixed_Genus-ASV_SILVA.biom"
                ),
            },
            "V1V3": {
                "input": (
                    "../Data/16S/Tables/"
                    "179426_feature-table_16S_V1V3_rare-11054_sampleIDfixed.biom"
                ),
                "output": (
                    "../Data/16S/Tables/"
                    "179426_feature-table_16S_V1V3_rare-11054_sampleIDfixed_Genus-ASV_SILVA.biom"
                ),
            },
        }

        for region, paths in biom_jobs.items():
            input_biom = paths["input"]
            output_biom = paths["output"]

            if not os.path.exists(input_biom):
                raise FileNotFoundError(input_biom)

            print(f"\nProcessing {region} table (SILVA)...")
            logging.info(f"Processing {region}: {input_biom}")

            df_counts = add_unique_tax_labels(
                input_biom,
                level="Genus"
            )

            save_biom_table(df_counts, output_biom)
            print(f"  → Saved {output_biom}")

        logging.info("All SILVA BIOM tables processed successfully.")
        print("\nDone.")

    except Exception as e:
        logging.error(f"Error occurred: {e}")
        raise
