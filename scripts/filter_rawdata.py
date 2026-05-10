import pandas as pd
from Bio import SeqIO
import argparse
import re
import numpy as np
import pathlib

segments = ["M", "L", "S"]
min_length_dic = {"S": 1000, "M": 3000, "L": 6000}



def rename_and_filter(df, exclude_accessions=None, allow_missing_date=False):
    # Only keep sequences with appropriate lengths
    df_filtered = df.dropna(subset=["Segment"])
    if not allow_missing_date:
        df_filtered = df_filtered.dropna(subset=["Isolate Collection date"])
    if exclude_accessions:
        df_filtered = df_filtered[~df_filtered["Accession"].isin(exclude_accessions)]
    df_filtered = df_filtered.loc[df_filtered["Length"] > 250]
    for segment in segments:
        min_length = min_length_dic[segment]
        df_filtered = df_filtered.drop(
            df_filtered[
                (df_filtered["Segment"] == segment)
                & (df_filtered["Length"] < min_length)
            ].index
        )
    # Add min_length_pass column
    df_filtered["min_length_pass"] = df_filtered.apply(
        lambda row: row["Length"] >= min_length_dic.get(row["Segment"], 0),
        axis=1,
    )
    # Rename columns
    df_filtered["country"] = df_filtered["Geographic Location"].apply(
        lambda x: x.split(":")[0] if isinstance(x, str) else None
    )
    df_renamed = df_filtered.rename(
        columns={
            "Virus Name": "virus",
            "Accession": "accession",
            "Isolate Collection date": "date",
            "Geographic Region": "region",
            "Submitter Names": "author"
        }
    )
    df_renamed["unique_id"] = df_renamed.apply(make_unique_name, axis=1)
    # Ensure uniqueness by adding -1, -2, ... if needed
    df_renamed = make_unique_names_within_segment(df_renamed, name_col="unique_id", segment_col="Segment")
    
    #counts = {}
    #unique_names = []
    #for name in df_renamed["unique_id"]:
    #    if name not in counts:
    #        counts[name] = 1
    #        unique_names.append(name)
    #    else:
    #        counts[name] += 1
    #        unique_names.append(f"{name}-{counts[name]}")
    #df_renamed["unique_id"] = unique_names

    return df_renamed

def group_name(isolate, date, country):
    return (
        str(isolate).replace(" ", "-")
        + "/"
        + str(date)
        + "/"
        + str(country).replace(" ", "-")
    )


def group_metadata(df):
    # Group sequences according to isolate and collection date
    df_grouped = df
    grouped = df_grouped.groupby(
        ["Isolate Lineage", "date", "country"]
    )
    groups = grouped.groups.keys()

    group_id = "None"
    df_grouped["group_id"] = group_id
    number_of_groups = 0
    for g in groups:
        isolate, date, location = g
        df_grouped["group_id"] = df_grouped.apply(
            lambda row: group_name(isolate, date, location)
            if row["Isolate Lineage"] == isolate
            and row["date"] == date
            and row["country"] == location
            else row["group_id"],
            axis=1,
        )
        number_of_groups += 1

    print("Number of groups: ", number_of_groups)

    # Remove groups with group_id = "None" -> no isolate given
    df_grouped = df_grouped.loc[df_grouped["group_id"] != "None"]

    # Add tag if all segments are present
    #all_segments = 0
    #for g in groups:
    #    isolate, date, location = g
    #    group_id = group_name(isolate, date, location)
    #    group_g = df_grouped.loc[df_grouped["group_id"] == group_id]
    #    if (
    #        "S" in group_g["Segment"].values
    #        and "M" in group_g["Segment"].values
    #        and "L" in group_g["Segment"].values
    #        and len(group_g) == 3
    #    ):
    #        df_grouped.loc[df_grouped["group_id"] == group_id, "nr_segments"] = "all"
    #        all_segments += 1
    #print("Number of groups with all segments: ", all_segments)

    #df_grouped = df_grouped.loc[df_grouped["nr_segments"] == "all"]

    # Write to directory: metadata only containing sequences where all segments are present
    path = pathlib.Path("data")
    path.mkdir(parents=True, exist_ok=True)
    df_grouped.to_csv("data/all_sequences_grouped.tsv", sep="\t", index=False)
    return df_grouped


def write_segment_metadata(df_all):
    for segment in segments:
        df_segment = df_all.loc[df_all["Segment"] == segment]
        df_segment.to_csv(
            "data/metadata_insdc_{0}.tsv".format(segment), sep="\t", index=False
        )


def write_fasta(all_sequences_path, segment, df_all, no_group=False):
    metadata = df_all.loc[df_all["Segment"] == segment]
    sequences_segment_path = "data/sequences_insdc_{0}.fasta".format(segment)
    with open(sequences_segment_path, "w") as seq:
        records = SeqIO.parse(all_sequences_path, "fasta")
        for record in records:
            corresponding_metadata = metadata.loc[metadata["accession"] == record.id]
            if len(corresponding_metadata) == 0:
                continue
            if no_group:
                # Use unique_id as sequence name
                unique_ids = list(corresponding_metadata["unique_id"])
                if len(unique_ids) == 0:
                    continue
                record.description = record.id
                record.id = unique_ids[0]
                SeqIO.write(record, seq, "fasta")
                continue
            group_ids = list(corresponding_metadata["group_id"])
            if len(group_ids) == 0:
                continue
            name = str(group_ids[0])
            record.description = record.id
            record.id = name
            SeqIO.write(record, seq, "fasta")


def write_segment_fasta(raw_sequences, df_all, no_group=False):
    all_sequences_path = "data/all_sequences_renamed.fasta"

    with open(all_sequences_path, "w") as renamed:
        records = SeqIO.parse(raw_sequences, "fasta")
        for record in records:
            record.description = record.id
            SeqIO.write(record, renamed, "fasta")
    for segment in segments:
        write_fasta(all_sequences_path, segment, df_all, no_group)

def make_unique_name(row):
    isolate = str(row.get("Isolate Lineage", "unknown")).replace(" ", "-")
    date = str(row.get("date", "unknown"))
    country = str(row.get("country", "unknown")).replace(" ", "-")
    return f"{isolate}/{date}/{country}"

# Ensure uniqueness by adding -1, -2, ... if needed, but only within each segment
def make_unique_names_within_segment(df, name_col="unique_id", segment_col="Segment"):
    df = df.copy()
    new_names = []
    for segment, group in df.groupby(segment_col):
        counts = {}
        for idx, name in group[name_col].items():
            if name not in counts:
                counts[name] = 1
                new_names.append((idx, name))
            else:
                counts[name] += 1
                new_names.append((idx, f"{name}-{counts[name]}"))
    # Assign new names back to the DataFrame
    name_map = dict(new_names)
    df[name_col] = df.index.map(name_map)
    return df

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="filter raw data",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "--metadata", type=str, required=True, help="csv file containing all sequences"
    )
    parser.add_argument(
        "--sequences",
        type=str,
        required=True,
        help="fasta file containing all sequences",
    )
    parser.add_argument(
        "--no-group",
        type=bool,
        default=False,
        help="whether to group sequences by isolate or not",
    )
    parser.add_argument(
        "--exclude-list",
        type=str,
        default=None,
        help="File with accession numbers to exclude (one per line)",
    )
    parser.add_argument(
        "--allow-missing-date",
        action="store_true",
        help="Allow sequences without Isolate Collection date",
    )

    args = parser.parse_args()
    exclude_accessions = set()
    if args.exclude_list:
        with open(args.exclude_list) as f:
            exclude_accessions = set(line.strip() for line in f if line.strip())

    df_metadata = pd.read_csv(args.metadata, sep="\t", on_bad_lines="warn")
    print("Loaded metadata:", "MN258226.1" in df_metadata["Accession"].values)
    df_renamed = rename_and_filter(df_metadata, exclude_accessions=exclude_accessions, allow_missing_date=args.allow_missing_date)
    print("After rename_and_filter:", "MN258226.1" in df_renamed["accession"].values)

    if args.no_group:
        write_segment_metadata(df_renamed)
        write_segment_fasta(args.sequences, df_renamed, args.no_group)
        exit()

    df_grouped = group_metadata(df_renamed)
    print("After group_metadata:", "MN258226.1" in df_grouped["accession"].values)

    write_segment_metadata(df_grouped)
    write_segment_fasta(args.sequences, df_grouped)
