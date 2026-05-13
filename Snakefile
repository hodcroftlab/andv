SEGMENTS = ["S", "M", "L"]

external_metadata = "" #"data/external_metadata.tsv"
external_fastas = {
    "S": "", #"data/external_fasta_S.fasta",
    "M": "", #"data/external_fasta_M.fasta",
    "L": "", #"data/external_fasta_L.fasta",
}


dropped_strains = ("config/dropped_strains.txt",)
reference = ("config/outgroup_{Segment}.gb",)
colors = ("config/colors.tsv",)
lat_longs = ("config/lat_longs.tsv",)
auspice_config = "config/auspice_config.json"
TAXON_ID = 1980456
LAPIS_URL = "https://lapis.pathoplexus.org/andv"

wildcard_constraints:
    Segment="S|M|L"

if os.uname().sysname == "Darwin":
    # Don't use conda-forge unzip on macOS
    # Due to https://github.com/conda-forge/unzip-feedstock/issues/16
    unzip = "/usr/bin/unzip"
else:
    unzip = "unzip"


# to clean input and force re-run, run as: 
# snakemake --cores 4 force_refresh && snakemake --cores 4 refresh_data && snakemake --cores 4 all

rule force_refresh:
    input:
        sequences=expand("data/sequences_pathoplexus_{Segment}.fasta", Segment=SEGMENTS),
        metadata=expand("data/metadata_pathoplexus_{Segment}.tsv", Segment=SEGMENTS),
    shell:
        """
        rm -f {input.sequences} \
              {input.metadata}
        """

rule refresh_data:
    input:
        sequences=expand("data/sequences_pathoplexus_{Segment}.fasta", Segment=SEGMENTS),
        metadata=expand("data/metadata_pathoplexus_{Segment}.tsv", Segment=SEGMENTS),

rule all:
    input:
        expand("auspice/andv_{Segment}.json", Segment=SEGMENTS),

rule fetch_pathoplexus_sequences:
    output:
        sequences="data/sequences_pathoplexus_{Segment}.fasta",
    params:
        lapis_url=LAPIS_URL,
    retries: 1
    shell:
        """
        curl -fsSL \
            '{params.lapis_url}/sample/unalignedNucleotideSequences/{wildcards.Segment}?versionStatus=LATEST_VERSION&isRevocation=false&downloadAsFile=true' \
            -o {output.sequences}
        """


rule fetch_pathoplexus_metadata:
    output:
        metadata="data/metadata_pathoplexus_{Segment}.tsv",
    params:
        lapis_url=LAPIS_URL,
    retries: 1
    shell:
        """
        curl -fsSL \
            '{params.lapis_url}/sample/details?versionStatus=LATEST_VERSION&isRevocation=false&dataFormat=TSV&downloadAsFile=true&length_{wildcards.Segment}From=1' \
            -o {output.metadata}
        """


rule add_external_sequences:
    input:
        ppx_metadata="data/metadata_pathoplexus_{Segment}.tsv",
        ppx_fasta="data/sequences_pathoplexus_{Segment}.fasta",
        external_metadata=lambda wildcards: external_metadata if external_metadata else [],
        external_fasta=lambda wildcards: external_fastas[wildcards.Segment] if external_fastas[wildcards.Segment] else [],
    output:
        metadata="data/metadata_{Segment}.tsv",
        fasta="data/sequences_{Segment}.fasta",
    run:
        import shutil
        import pandas as pd
        # Concatenate metadata
        ncbi_meta = pd.read_csv(input.ppx_metadata, sep="\t")
        if input.external_metadata:
            ext_meta = pd.read_csv(input.external_metadata, sep="\t")
            # Copy external metadata for each sequence in external fasta
            from Bio import SeqIO
            ext_fasta_records = list(SeqIO.parse(input.external_fasta, "fasta")) if input.external_fasta else []
            ext_meta_rows = []
            for record in ext_fasta_records:
                row = ext_meta.iloc[0].copy()
            #    row["accession"] = record.id
                ext_meta_rows.append(row)
            ext_meta_expanded = pd.DataFrame(ext_meta_rows)
            combined_meta = pd.concat([ncbi_meta, ext_meta_expanded], ignore_index=True)
        else:
            combined_meta = ncbi_meta
        combined_meta.to_csv(output.metadata, sep="\t", index=False)
        # Concatenate fasta
        with open(output.fasta, "wb") as out_f:
            with open(input.ppx_fasta, "rb") as in1:
                shutil.copyfileobj(in1, out_f)
            if input.external_fasta:
                with open(input.external_fasta, "rb") as in2:
                    out_f.write(b"\n")
                    shutil.copyfileobj(in2, out_f)

rule format_metadata:
    message:
        "Formatting metadata for the run"
    input:
        metadata="data/metadata_{Segment}.tsv",
    output:
        metadata="data/metadata_formatted_{Segment}.tsv",
    shell:
        """
        python scripts/format_metadata.py \
            --metadata {input.metadata} \
            --output {output.metadata} \
            --segment {wildcards.Segment}
        """


rule index_sequences:
    message:
        """
        Creating an index of sequence composition for filtering.
        """
    input:
        sequences="data/sequences_{Segment}.fasta",
    output:
        sequence_index="results/sequence_index_{Segment}.tsv",
    shell:
        """
        augur index \
            --sequences {input.sequences} \
            --output {output.sequence_index}
        """


rule filter:
    message:
        """
        Filtering to
          - {params.sequences_per_group} sequence(s) per {params.group_by!s}
          - from {params.min_date} onwards
          - excluding strains in {input.exclude}
        """
    input:
        sequences="data/sequences_{Segment}.fasta",
        sequence_index="results/sequence_index_{Segment}.tsv",
        metadata="data/metadata_formatted_{Segment}.tsv",
        exclude=dropped_strains,
    output:
        filtered_metadata="results/initial_filtered_metadata_{Segment}.tsv",
        filtered_sequences="results/initial_filtered_sequences_{Segment}.fasta",
    params:
        group_by="country year",
        sequences_per_group=20,
        id_column="accessionVersion",
        min_date=1900,
        min_length=lambda wildcards: {"S": 1000, "M": 3000, "L": 6000}[wildcards.Segment],
        #query="min_length_pass == 'True'",
        #query="nr_segments == 'all'",
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --sequence-index {input.sequence_index} \
            --metadata {input.metadata} \
            --exclude {input.exclude} \
            --output-metadata {output.filtered_metadata} \
            --metadata-id-columns {params.id_column} \
            --min-length {params.min_length} \
            --output-sequences {output.filtered_sequences}
        """
        #             --group-by {params.group_by} \
        #    --sequences-per-group {params.sequences_per_group} \
        #          --min-date {params.min_date} 
        #            --query {params.query:q}



rule align:
    message:
        """
        Aligning sequences to {input.reference}
          - filling gaps with N
        """
    input:
        sequences=rules.filter.output.filtered_sequences,
        reference=reference,
    output:
        alignment="results/aligned_{Segment}.fasta",
    shell:
        """
        augur align \
            --sequences {input.sequences} \
            --reference-sequence {input.reference} \
            --remove-reference\
            --output {output.alignment} \
            --fill-gaps
        """


rule tree:
    message:
        "Building tree"
    input:
        alignment=rules.align.output.alignment,
    output:
        tree="results/tree_raw_{Segment}.nwk",
    params:
        tree_builder_args="-czb -st DNA",
    shell:
        """
        augur tree \
            --alignment {input.alignment} \
            --tree-builder-args="{params.tree_builder_args}" \
            --output {output.tree}
        """
        # Added -st DNA because it was unable to detect sequence type automatically
        # -czb is a flag that collapses zero branches 


rule refine:
    message:
        """
        Refining tree
          - estimate timetree
          - use {params.coalescent} coalescent timescale
          - estimate {params.date_inference} node dates
          - filter tips more than {params.clock_filter_iqd} IQDs from clock expectation
        """
    input:
        tree=rules.tree.output.tree,
        alignment=rules.align.output,
        metadata="data/metadata_formatted_{Segment}.tsv",
    output:
        tree="results/tree_{Segment}.nwk",
        node_data="results/branch_lengths_{Segment}.json",
    params:
        coalescent="opt",
        date_inference="marginal",
        clock_filter_iqd=4,
        id_column="unique_id",
        root="mid_point",
    shell:
        """
        augur refine \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --metadata {input.metadata} \
            --output-tree {output.tree} \
            --output-node-data {output.node_data} \
            --root {params.root} \
            --metadata-id-columns {params.id_column}
        """

rule ancestral:
    message:
        "Reconstructing ancestral sequences and mutations"
    input:
        tree="results/tree_{Segment}.nwk",
        alignment=rules.align.output,
        reference=reference,
    output:
        node_data="results/nt_muts_{Segment}.json",
    params:
        inference="joint",
    shell:
        """
        augur ancestral \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --output-node-data {output.node_data} \
            --root-sequence {input.reference} \
            --inference {params.inference}
        """


rule translate:
    message:
        "Translating amino acid sequences"
    input:
        tree="results/tree_{Segment}.nwk",
        node_data=rules.ancestral.output.node_data,
        reference=reference,
    output:
        node_data="results/aa_muts_{Segment}.json",
    shell:
        """
        augur translate \
            --tree {input.tree} \
            --ancestral-sequences {input.node_data} \
            --reference-sequence {input.reference} \
            --output-node-data {output.node_data} \
        """


rule clades:
    input:
        tree=rules.refine.output.tree,
    output:
        node_data="results/clades_{Segment}.json",
    shell:
        """
        python scripts/get_clades.py \
            --tree {input.tree} \
            --node-data {output.node_data} \
            --clade-name {wildcards.Segment}
        """


rule traits:
    message:
        "Inferring ancestral traits for {params.columns!s}"
    input:
        tree="results/tree_{Segment}.nwk",
        metadata="data/metadata_formatted_{Segment}.tsv",
    output:
        node_data="results/traits_{Segment}.json",
    params:
        columns="region country",
        id_column="group_id",
    shell:
        """
        augur traits \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --output-node-data {output.node_data} \
            --columns {params.columns} \
            --confidence \
            --metadata-id-columns {params.id_column}
        """


rule export:
    message:
        "Exporting data files for for auspice"
    input:
        tree="results/tree_{Segment}.nwk",
        metadata="data/metadata_formatted_{Segment}.tsv",
        clades=rules.clades.output.node_data,
        branch_lengths=rules.refine.output.node_data,
        #traits=rules.traits.output.node_data,
        nt_muts=rules.ancestral.output.node_data,
        aa_muts=rules.translate.output.node_data,
        colors=colors,
        lat_longs=lat_longs,
        auspice_config="config/auspice_config_{Segment}.json",
        description = "config/description.md"
    output:
        auspice_json="auspice/andv_{Segment}.json",
    params:
        id_column="accessionVersion",
        metadata_columns="unique_id dataUseTerms restrictedUntil PPX_accession INSDC_accession" #dataUseTerms__url
    shell:
        """
        augur export v2 \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --metadata-columns {params.metadata_columns} \
            --node-data {input.branch_lengths} {input.nt_muts} {input.aa_muts} {input.clades} \
            --colors {input.colors} \
            --lat-longs {input.lat_longs} \
            --auspice-config {input.auspice_config} \
            --output {output.auspice_json} \
            --metadata-id-columns {params.id_column} \
            --include-root-sequence-inline \
            --description "{input.description}" \
            --metadata-id-columns "unique_id"
        """


rule clean:
    message:
        "Removing directories: {params}"
    params:
        "results ",
        "auspice",
    shell:
        "rm -rfv {params}"
