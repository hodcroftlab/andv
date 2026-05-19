include: "Snakefile"


rule genome:
    input:
        "auspice/andv_hondius-outbreak-genome.json"


# Looking across all 3 segment trees this clade (c. 20 strains) doesn't appear to have any reassortment
# (at least, none that's disrupted the tree). 
NON_REASSORTING_CLADE = [
    'PP_006X9ND.1', # nan/2026-05-10/France-1 (Hondius outbreak)
    'PP_006W2KU.1', # CHI-Hu13724_P2/2024-07-30/Chile
]

rule outbreak_samples:
    """For each segment find the tips previously classified as part of the outbreak via `augur clades`"""
    input:
        tree=expand("results/tree_{Segment}.nwk", Segment=SEGMENTS),
    params:
        segments=SEGMENTS,
        ca=NON_REASSORTING_CLADE,
    output:
        include="results/genome/samples_genome.txt",
    run:
        from Bio import Phylo
        segment_samples = []
        for segment, fname in zip(params.segments, input.tree):
            tree = Phylo.read(fname, "newick")
            mrca = tree.common_ancestor({"name": params.ca[0]}, {"name": params.ca[1]})
            samples = {tip.name for tip in mrca.get_terminals()}
            print(f"Segment {segment} n={len(samples)} samples")
            segment_samples.append(samples)
        common_samples = set.intersection(*segment_samples)
        print(f"Samples in common: n={len(common_samples)} samples")
        with open(output.include, 'w') as fh:
            for sample in sorted(common_samples):
                print(sample, file=fh)

rule alignment_genome_segment:
    """Pull out the segment alignments for the outbreak samples"""
    input:
        alignment="results/aligned_{Segment}.fasta",
        include="results/genome/samples_genome.txt",
        metadata="data/metadata_curated_{Segment}.tsv"
    output:
        alignment="results/genome/aligned_{Segment}.fasta",
        metadata="results/genome/metadata_{Segment}.tsv"
    shell:
        r"""
        augur filter \
            --metadata-id-columns accessionVersion \
            --sequences {input.alignment} --metadata {input.metadata} \
            --exclude-all --include {input.include} \
            --output-sequences {output.alignment} --output-metadata {output.metadata}
        """
    
rule genome_alignment:
    """Concatenate segment alignments"""
    input:
        alignments=expand("results/genome/aligned_{Segment}.fasta", Segment=SEGMENTS), # order: canonical S,M,L
    output:
        alignment="results/genome/aligned.fasta"
    shell:
        r"""
        seqkit concat {input.alignments} > {output.alignment}
        """

rule genome_masking:
    """Mask terminal regions of segments"""
    input:
        alignment="results/genome/aligned.fasta",
        mask="config/genome-mask.bed",
    output:
        masked="results/genome/masked.fasta",
    shell:
        r"""
        augur mask --sequences {input.alignment} --mask {input.mask} --output {output.masked}
        """

rule genome_tree:
    input:
        alignment="results/genome/masked.fasta",
    output:
        tree="results/genome/tree_raw.nwk",
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


rule genome_refine:
    input:
        tree="results/genome/tree_raw.nwk",
        alignment="results/genome/masked.fasta",
        metadata="results/genome/metadata_S.tsv" # S segment picked randomly
    output:
        tree="results/genome/tree.nwk",
        node_data="results/genome/branch_lengths.json",
    params:
        id_column="accessionVersion",
        root = "mid_point",
    shell:
        """
        augur refine \
            --metadata-id-columns {params.id_column} \
            --metadata {input.metadata} \
            --alignment {input.alignment} \
            --tree {input.tree} \
            --output-tree {output.tree} \
            --output-node-data {output.node_data} \
            --root {params.root} \
        """

rule ancestral_genome:
    input:
        tree="results/genome/tree.nwk",
        alignment="results/genome/aligned.fasta", # not masked, so we'll report mutations in masked regions
    output:
        node_data="results/genome/nt_muts.json",
    params:
        inference="joint",
    shell:
        """
        augur ancestral \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --output-node-data {output.node_data} \
            --inference {params.inference}
        """


rule translate_genome:
    input:
        tree="results/genome/tree.nwk",
        node_data="results/genome/nt_muts.json",
        reference="config/genome.gb",
    output:
        node_data="results/genome/aa_muts.json",
    shell:
        """
        augur translate \
            --tree {input.tree} \
            --ancestral-sequences {input.node_data} \
            --reference-sequence {input.reference} \
            --output-node-data {output.node_data} \
        """


rule clades_genome:
    input:
        tree="results/genome/tree.nwk",
        nuc_muts = "results/genome/nt_muts.json",
        aa_muts = "results/genome/aa_muts.json",
        clades = "config/clades_S.tsv", # numbering good as S is starting segment
    output:
        node_data="results/genome/clades.json",
    shell:
        """
        augur clades --tree {input.tree} \
            --mutations {input.nuc_muts} {input.aa_muts} \
            --clades {input.clades} \
            --output-node-data {output.node_data}
        """

rule traits_genome:
    input:
        tree="results/genome/tree.nwk",
        metadata="results/genome/metadata_S.tsv", # S segment picked randomly
    output:
        node_data="results/genome/traits.json",
    params:
        columns=["country"],
        id_column="accessionVersion",
    shell:
        """
        augur traits \
            --tree {input.tree} \
            --metadata-id-columns {params.id_column} \
            --metadata {input.metadata} \
            --output-node-data {output.node_data} \
            --columns {params.columns} \
            --confidence
        """


rule export_genome:
    input:
        auspice_config="config/auspice_config_genome.json",
        tree="results/genome/tree.nwk",
        metadata="results/genome/metadata_S.tsv", # S segment picked randomly
        node_data=[
            "results/genome/branch_lengths.json",
            # "results/genome/traits.json", # skip country DTA as it's misleading for the hondius outbreak
            "results/genome/nt_muts.json",
            "results/genome/aa_muts.json",
            "results/genome/clades.json",
        ],
        lat_longs=lat_longs,
        description = "config/description.md"
    output:
        auspice_json="auspice/andv_hondius-outbreak-genome.json",
    params:
        id_column="accessionVersion",
        metadata_columns="unique_id dataUseTerms restrictedUntil PPX_accession INSDC_accession hostNameCommon hostNameScientific" #dataUseTerms__url
    shell:
        """
        augur export v2 \
            --metadata-id-columns {params.id_column} \
            --metadata {input.metadata} \
            --tree {input.tree} \
            --auspice-config {input.auspice_config} \
            --lat-longs {input.lat_longs} \
            --metadata-columns {params.metadata_columns} \
            --node-data {input.node_data} \
            --description {input.description} \
            --include-root-sequence-inline \
            --output {output.auspice_json}
        """
