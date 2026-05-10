SEGMENTS = ["M", "L", "S"]

external_metadata = "data/external_metadata.tsv"
external_fastas = {
    "S": "data/external_fasta_S.fasta",
    "M": "data/external_fasta_M.fasta",
    "L": "data/external_fasta_L.fasta",
}

rule all:
    input:
        expand("auspice/andv_{Segment}.json", Segment=SEGMENTS),


dropped_strains = ("config/dropped_strains.txt",)
reference = ("config/outgroup_{Segment}.gb",)
colors = ("config/colors.tsv",)
lat_longs = ("config/lat_longs.tsv",)
auspice_config = "config/auspice_config.json"
TAXON_ID = 1980456

if os.uname().sysname == "Darwin":
    # Don't use conda-forge unzip on macOS
    # Due to https://github.com/conda-forge/unzip-feedstock/issues/16
    unzip = "/usr/bin/unzip"
else:
    unzip = "unzip"

rule fetch_ncbi_dataset_package:
    output:
        dataset_package="data/ncbi_dataset.zip",
    shell:
        """
        datasets download virus genome taxon {TAXON_ID} \
            --no-progressbar \
            --filename {output.dataset_package} \
        """


rule extract_ncbi_dataset_sequences:
    input:
        dataset_package="data/ncbi_dataset.zip",
    output:
        ncbi_dataset_sequences="data/sequences.fasta",
    params:
        unzip=unzip,
    shell:
        """
        {params.unzip} -o {input.dataset_package} -d data 
        seqkit seq -w0 data/ncbi_dataset/data/genomic.fna \
        > {output.ncbi_dataset_sequences}
        """


rule format_ncbi_dataset_report:
    input:
        dataset_package="data/ncbi_dataset.zip",
    output:
        ncbi_dataset_tsv="data/metadata_post_extract.tsv",
    shell:
        """
        dataformat tsv virus-genome \
            --package {input.dataset_package} \
            > {output.ncbi_dataset_tsv} 
        """

checkpoint filter_rawdata:
    input:
        sequences="data/sequences.fasta",
        ncbi_dataset_tsv="data/metadata_post_extract.tsv",
    output:
        sequences=expand("data/sequences_insdc_{Segment}.fasta", Segment=SEGMENTS),
        metadata=expand("data/metadata_insdc_{Segment}.tsv", Segment=SEGMENTS),
    shell:
        """
        python scripts/filter_rawdata.py \
            --sequences {input.sequences} \
            --metadata {input.ncbi_dataset_tsv} \
            --exclude-list "config/rawdataexclude.txt" \
            --no-group true --allow-missing-date
        """

rule add_external_sequences:
    input:
        ncbi_metadata="data/metadata_insdc_{Segment}.tsv",
        ncbi_fasta="data/sequences_insdc_{Segment}.fasta",
        external_metadata=lambda wildcards: external_metadata if external_metadata else [],
        external_fasta=lambda wildcards: external_fastas[wildcards.Segment] if external_fastas[wildcards.Segment] else [],
    output:
        metadata="data/metadata_{Segment}.tsv",
        fasta="data/sequences_{Segment}.fasta",
    run:
        import shutil
        import pandas as pd
        # Concatenate metadata
        ncbi_meta = pd.read_csv(input.ncbi_metadata, sep="\t")
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
            with open(input.ncbi_fasta, "rb") as in1:
                shutil.copyfileobj(in1, out_f)
            if input.external_fasta:
                with open(input.external_fasta, "rb") as in2:
                    out_f.write(b"\n")
                    shutil.copyfileobj(in2, out_f)

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
        metadata="data/metadata_{Segment}.tsv",
        exclude=dropped_strains,
    output:
        filtered_metadata="results/initial_filtered_metadata_{Segment}.tsv",
    params:
        group_by="country year",
        sequences_per_group=20,
        id_column="unique_id",
        min_date=1900,
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
            --metadata-id-columns {params.id_column} 
        
        """
        #             --group-by {params.group_by} \
        #    --sequences-per-group {params.sequences_per_group} \
        #          --min-date {params.min_date} 
        #            --query {params.query:q}


rule extract_selected_groups:
    message:
        """
        Select the group id of filtered sequences
        """
    input:
        metadata=rules.filter.output.filtered_metadata,
    output:
        filtered_groups_list="results/filtered_groups_{Segment}.tsv",
        filtered_groups="results/filtered_groups_{Segment}.txt",
    shell:
        """
        tsv-select -H -f unique_id {input.metadata} > {output.filtered_groups_list}
        cp {output.filtered_groups_list} {output.filtered_groups}
        """


rule select_sequences:
    message:
        """
        Select all sequences according to group-id
        """
    input:
        groups=rules.extract_selected_groups.output.filtered_groups,
        sequences="data/sequences_{Segment}.fasta",
    output:
        selected_sequences="results/selected_sequences_{Segment}.fasta",
    shell:
        """
        seqkit grep -f {input.groups} < {input.sequences} > {output.selected_sequences}
        """


rule align:
    message:
        """
        Aligning sequences to {input.reference}
          - filling gaps with N
        """
    input:
        sequences=rules.select_sequences.output.selected_sequences,
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
    shell:
        """
        augur tree \
            --alignment {input.alignment} \
            --tree-builder-args="-czb" \
            --output {output.tree}
        """


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
        metadata="data/metadata_{Segment}.tsv",
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


rule use_treeknit:
    message:
        """
        Using treeknit to detangle trees
        """
    input:
        tree_S="results/tree_S.nwk",
        tree_M="results/tree_M.nwk",
        tree_L="results/tree_L.nwk",
    output:
        tree_S_res="results/tree_S_resolved.nwk",
        tree_M_res="results/tree_M_resolved.nwk",
        tree_L_res="results/tree_L_resolved.nwk",
    shell:
        """
        treeknit {input.tree_S} {input.tree_L} {input.tree_M} --better-trees --rounds 3 -o results
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
        metadata="data/metadata_{Segment}.tsv",
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
        metadata="data/metadata_{Segment}.tsv",
        clades=rules.clades.output.node_data,
        branch_lengths=rules.refine.output.node_data,
        #traits=rules.traits.output.node_data,
        nt_muts=rules.ancestral.output.node_data,
        aa_muts=rules.translate.output.node_data,
        colors=colors,
        lat_longs=lat_longs,
        auspice_config="config/auspice_config_{Segment}.json",
    output:
        auspice_json="auspice/andv_{Segment}.json",
    params:
        id_column="unique_id",
    shell:
        """
        augur export v2 \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --metadata-columns accession \
            --node-data {input.branch_lengths} {input.nt_muts} {input.aa_muts} {input.clades} \
            --colors {input.colors} \
            --lat-longs {input.lat_longs} \
            --auspice-config {input.auspice_config} \
            --output {output.auspice_json} \
            --metadata-id-columns {params.id_column}
        """


rule clean:
    message:
        "Removing directories: {params}"
    params:
        "results ",
        "auspice",
    shell:
        "rm -rfv {params}"
