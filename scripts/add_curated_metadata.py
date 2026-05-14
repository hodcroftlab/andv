#!/usr/bin/env python3
"""
Add curated metadata from ARTIC ANDV repository to main metadata.

This script merges curated metadata (from andv-metadata.csv) with the main
metadata table, matching by strain/accession identifiers and segment.
"""

import argparse
import pandas as pd


def main():
    parser = argparse.ArgumentParser(
        description="Add curated metadata to main metadata table"
    )
    parser.add_argument(
        "--main-metadata",
        required=True,
        help="Path to main metadata TSV file"
    )
    parser.add_argument(
        "--curated-metadata",
        required=True,
        help="Path to curated metadata CSV file from ARTIC ANDV"
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Path to output metadata TSV file"
    )
    parser.add_argument(
        "--segment",
        required=True,
        choices=["S", "M", "L"],
        help="Segment being processed (S, M, or L)"
    )

    args = parser.parse_args()

    # Read input files
    main_meta = pd.read_csv(args.main_metadata, sep="\t")
    curated_meta = pd.read_csv(args.curated_metadata)

    # Map segment to accession column names
    # Main metadata has columns: insdcAccessionFull_S, insdcAccessionFull_M, insdcAccessionFull_L
    accession_column = f"insdcAccessionFull_{args.segment}"
    # Curated metadata has columns named: S, M, L (for segment-specific accessions)
    curated_accession_column = args.segment

    # Extract base accession (without version) from the versioned accession
    def get_base_accession(versioned_acc):
        if pd.isna(versioned_acc):
            return None
        # Remove version suffix (e.g., "AF291702.1" -> "AF291702")
        return str(versioned_acc).split(".")[0]

    # Create lookup dictionary from curated metadata
    curated_lookup = {}
    for idx, row in curated_meta.iterrows():
        acc = row[curated_accession_column]
        if pd.notna(acc) and acc != "":
            base_acc = get_base_accession(acc)
            curated_lookup[base_acc] = row

    # Let's an an empty column to main_meta to hold 'province' information 
    if "province" not in main_meta.columns:
        main_meta["province"] = pd.NA

    # First let's do country and region (province)- but only if the main metadata doesn't already have this information (i.e. if it's NA)
    for idx, row in main_meta.iterrows():
        versioned_acc = row[accession_column]
        if pd.notna(versioned_acc):
            base_acc = get_base_accession(versioned_acc)
            if base_acc in curated_lookup:
                curated_row = curated_lookup[base_acc]
                # only replace in main_meta if each is empty/NA in main_meta
                if pd.isna(row.get("country")) and pd.notna(curated_row["country"]):
                    main_meta.at[idx, "country"] = curated_row["country"]
                if pd.isna(row.get("province")) and pd.notna(curated_row["region"]):
                    # replace any _ with ' ', e.g. "San_Pedro" -> "San Pedro"
                    curated_region = str(curated_row["region"]).replace("_", " ")
                    # if curated_region is "ex chile", set to "" (this seems to refer to travel history, not more exact location)
                    if curated_region.lower() == "ex chile":
                        curated_region = ""
                    main_meta.at[idx, "province"] = curated_region

    # Now let's do date - only if date is empty in main_meta, get it from curated_meta
    for idx, row in main_meta.iterrows():
        versioned_acc = row[accession_column]
        if pd.notna(versioned_acc):
            base_acc = get_base_accession(versioned_acc)
            if base_acc in curated_lookup:
                curated_row = curated_lookup[base_acc]
                if pd.isna(row.get("date")) and pd.notna(curated_row["date"]):
                    # first, curate the curated_data date into the right format:
                    # if 1997 - should be 1997-XX-XX, if 1997-05 - should be 1997-05-XX, if 1997-05-12 - should be 1997-05-12
                    curated_date = str(curated_row["date"])
                    if len(curated_date) == 4:  # year only
                        curated_date = f"{curated_date}-XX-XX"
                    elif len(curated_date) == 7:  # year and month
                        curated_date = f"{curated_date}-XX"
                    main_meta.at[idx, "date"] = curated_date
                else:
                    #if not empty, double-check it's in the right format, 
                    # if 1997 - should be 1997-XX-XX, if 1997-05 - should be 1997-05-XX, if 1997-05-12 - should be 1997-05-12
                    existing_date = str(row.get("date"))
                    if len(existing_date) == 4:  # year only
                        main_meta.at[idx, "date"] = f"{existing_date}-XX-XX"
                    elif len(existing_date) == 7:  # year and month
                        main_meta.at[idx, "date"] = f"{existing_date}-XX"   

    # Finally do the host - only if hostNameScientific is empty in main_meta, get it from curated_meta
    for idx, row in main_meta.iterrows():
        versioned_acc = row[accession_column]
        if pd.notna(versioned_acc):
            base_acc = get_base_accession(versioned_acc)
            if base_acc in curated_lookup:
                curated_row = curated_lookup[base_acc]
                if pd.isna(row.get("hostNameScientific")) and pd.notna(curated_row["host"]):
                    # replace any _ in the name with a space, e.g. Oligoryzomys_longicaudatus -> Oligoryzomys longicaudatus
                    curated_host = str(curated_row["host"]).replace("_", " ")
                    main_meta.at[idx, "hostNameScientific"] = curated_host

    # Sometimes there is more info in the main_meta geoLocAdmin1 column
    # If that column is not empty, and province is empty, move the info from geoLocAdmin1 to province
    # But, if it has a comma ("Castelo dos Sonhos, Mato Grasso State")
    # separate by comma and take the second value (province rather than city)
    for idx, row in main_meta.iterrows():
        if pd.notna(row.get("geoLocAdmin1")) and pd.isna(row.get("province")):
            geo_loc = str(row["geoLocAdmin1"])
            if "," in geo_loc:
                main_meta.at[idx, "province"] = geo_loc.split(",")[1].strip()
            else:
                main_meta.at[idx, "province"] = geo_loc

    # Do some semi-manual fixing on some of the host columns.
    # If hostNameScientific is "Homo sapiens" and hostNameCommon is empty, set hostNameCommon to "human"
    for idx, row in main_meta.iterrows():
        if row.get("hostNameScientific") == "Homo sapiens" and pd.isna(row.get("hostNameCommon")):
            main_meta.at[idx, "hostNameCommon"] = "human"

    ##### Manual fixes
    # define a list of manual fixes for hostNameCommon given hostNameScientific
    manual_fixes = {
        "A.azarae": ("Akodon azarae", "Azara's grass mouse"),
        "Oligoryzomys longicaudatus": ("Oligoryzomys longicaudatus", "long-tailed pygmy rice rat"),
        "Oxymycterus nasutus": ("Oxymycterus nasutus", "long-nosed hocicudo"),
        "Oligoryzomys flavescens": ("Oligoryzomys flavescens", "yellow pygmy rice rat"),
        "Akodon montensis": ("Akodon montensis", "montane grass mouse"),
        "Oligoryzomys nigripes": ("Oligoryzomys nigripes", "black-footed pygmy rice rat"),
        "Mesocricetus auratus": ("Mesocricetus auratus", "Syrian golden hamster"),
        "Oligoryzomyz flavescens": ("Oligoryzomys flavescens", "yellow pygmy rice rat"),
    }

    # make the fixes
    for idx, row in main_meta.iterrows():
        if row.get("hostNameScientific") in manual_fixes:
            main_meta.at[idx, "hostNameScientific"], main_meta.at[idx, "hostNameCommon"] = manual_fixes[row.get("hostNameScientific")]

    result = main_meta

    # Write output
    result.to_csv(args.output, sep="\t", index=False)
    print(f"Wrote curated metadata to {args.output}")


if __name__ == "__main__":
    main()
