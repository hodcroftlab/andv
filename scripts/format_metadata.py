import pandas as pd
from Bio import SeqIO
import argparse
import re
import numpy as np
import pathlib

#Make a unique, useful name that includes strain, date, and country
def make_unique_name(row):
    isolate = str(row.get("isolate", "unknown")).replace(" ", "-")
    date = str(row.get("date", "unknown"))
    country = str(row.get("country", "unknown")).replace(" ", "-")
    return f"{isolate}/{date}/{country}"

# ensure that names are not duplicated, if they are, add a number to the end of the name, e.g. name-1, name-2, etc.
def deduplicate_names(names):
    counts = {}
    for name in names:
        counts[name] = counts.get(name, 0) + 1
    occurrence = {}
    result = []
    for name in names:
        if counts[name] > 1:
            occurrence[name] = occurrence.get(name, 0) + 1
            result.append(f"{name}-{occurrence[name]}")
        else:
            result.append(name)
    return result

# Rename the following columns as such:
# sampleCollectionDate -> date
# geoLocCountry -> country
# specimenCollectorSampleId -> isolate
# accession -> PPX_accession
# dataUseTermsUrl -> dataUseTerms__url
# dataUseTermsRestrictedUntil -> restrictedUntil

def rename_columns(df):
	df = df.rename(
		columns={
			"sampleCollectionDate": "date",
			"geoLocCountry": "country",
			"specimenCollectorSampleId": "isolate",
			"accession": "PPX_accession",
			"dataUseTermsRestrictedUntil": "restrictedUntil",
            "dataUseTermsUrl": "dataUseTerms__url"
		}
	)
	return df

# For the INSDC ID we need to treat this according to segment - so pass in segment
#insdcAccessionBase_{Segment} -> INSDC_accession
def add_insdc_accession(df, segment):
	insdc_accession_col = f"insdcAccessionBase_{segment}"
	if insdc_accession_col in df.columns:
		df = df.rename(columns={insdc_accession_col: "INSDC_accession"})
	else:
		df["INSDC_accession"] = np.nan
	return df


#create URL fields - need to attach the accessions to URLs to make URL fields
# do something like this
#for index, record in enumerate(records):
#	record = record.copy()
#
#	ppx_accession = record.get('PPX_accession', None) 
#	insdc_accession = record.get('INSDC_accession', None) 
# Add INSDC_accession__url and PPX_accession__url fields to NDJSON records
#record['PPX_accession__url'] = f"https://pathoplexus.org/seq/{ppx_accession}" \
#	if ppx_accession \
#	else ""
#record['INSDC_accession__url'] = f"https://www.ncbi.nlm.nih.gov/nuccore/{insdc_accession}" \
#	if insdc_accession \
#	else ""
def add_url_fields(df):
    df = df.copy()
    df['PPX_accession__url'] = df['PPX_accession'].apply(lambda x: f"https://pathoplexus.org/seq/{x}" if x else "")
    df['INSDC_accession__url'] = df['INSDC_accession'].apply(lambda x: f"https://www.ncbi.nlm.nih.gov/nuccore/{x}" if x else "")
    return df

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="format metadata",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "--metadata", type=str, required=True, help="csv file containing all sequences"
    )
    parser.add_argument(
        "--output", type=str, required=True, help="csv file to write formatted metadata"
    )
    parser.add_argument(
        "--segment", type=str, required=True, help="segment being processed"
    )

    args = parser.parse_args()

    df = pd.read_csv(args.metadata, sep="\t")
    df = rename_columns(df)
    df["unique_id"] = deduplicate_names(df.apply(make_unique_name, axis=1))
    df = add_insdc_accession(df, args.segment)
    df = add_url_fields(df)
    
	# Write the formatted metadata to a new CSV file
    df.to_csv(args.output, sep="\t", index=False)
