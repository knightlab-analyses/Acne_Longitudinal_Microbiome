#!/usr/bin/env python

import os
import logging
import warnings
import biom
import pandas as pd
from biom import load_table
from biom.util import biom_open

######################################################
# COLLAPSE SILVA GENUS-ASV BIOM TABLES TO GENUS LEVEL
######################################################

warnings.simplefilter(action="ignore", category=FutureWarning)

# Logging
os.makedirs("../Logs", exist_ok=True)
logging.basicConfig(
    filename="../Logs/Genus-level_tbls.log",
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s"
)

def collapse_genus_asv_biom(biom_path: str) -> pd.DataFrame:
    """
    Collapse a SILVA Genus-ASV BIOM table to pure Genus level.

    Input:
        samples × Genus ASV (e.g. g__Cutibacterium ASV-1)

    Output:
        samples × Genus
    """

    logging.info(f"Loading BIOM: {biom_path}")
    table = load_table(biom_path)

    # samples × features
    df = table.to_dataframe().T
    logging.info(f"Input shape: {df.shape}")

    # Extract genus name by removing ASV suffix
    df.columns = (
        df.columns
        .str.replace(r"\s+ASV-\d+$", "", regex=True)
        .str.strip()
    )

    # Collapse by genus
    df_genus = df.groupby(df.columns, axis=1).sum()

    logging.info(f"Collapsed shape: {df_genus.shape}")
    return df_genus


def save_biom(df: pd.DataFrame, out_path: str):

    os.makedirs(os.path.dirname(out_path), exist_ok=True)

    table = biom.table.Table(
        df.T.values,
        observation_ids=df.columns,
        sample_ids=df.index
    )

    with biom_open(out_path, "w") as f:
        table.to_hdf5(
            f,
            generated_by="SILVA Genus collapse"
        )

    logging.info(f"Saved BIOM: {out_path}")


# Main
if __name__ == "__main__":

    biom_jobs = {
        "V4": {
            "input": (
                "../Data/16S/Tables/"
                "174951_feature-table_16S_V4_rare-3769_sampleIDfixed_Genus-ASV_SILVA.biom"
            ),
            "output": (
                "../Data/16S/Tables/"
                "174951_feature-table_16S_V4_rare-3769_sampleIDfixed_Genus_SILVA.biom"
            ),
        },
        "V1V3": {
            "input": (
                "../Data/16S/Tables/"
                "179426_feature-table_16S_V1V3_rare-11054_sampleIDfixed_Genus-ASV_SILVA.biom"
            ),
            "output": (
                "../Data/16S/Tables/"
                "179426_feature-table_16S_V1V3_rare-11054_sampleIDfixed_Genus_SILVA.biom"
            ),
        },
    }

    try:
        for region, paths in biom_jobs.items():

            print(f"\nCollapsing {region} SILVA table → Genus")
            logging.info(f"Processing {region}")

            if not os.path.exists(paths["input"]):
                raise FileNotFoundError(paths["input"])

            df_genus = collapse_genus_asv_biom(paths["input"])
            save_biom(df_genus, paths["output"])

            print(f"  → Saved {paths['output']}")

        print("\nDone.")
        logging.info("All SILVA Genus tables collapsed successfully.")

    except Exception as e:
        logging.error(f"Error: {e}")
        raise
