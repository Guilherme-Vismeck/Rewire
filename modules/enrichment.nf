process PREPARE_GENE_LISTS {
    tag "gene_lists"
    publishDir "${params.outdir}/enrichment", mode: 'copy', pattern: '*.{tsv,txt}'

    input:
    path de_vs_rewiring
    path stable
    path gene_modules

    output:
    path "gene_lists.tsv", emit: gene_lists
    path "universe.txt", emit: universe
    path "gene_list_sizes.tsv", emit: sizes

    script:
    """
    gene_lists.py \\
        --de_vs_rewiring ${de_vs_rewiring} \\
        --stable ${stable} \\
        --gene_modules ${gene_modules}
    """
}

process ENRICH_CUSTOM {
    tag "custom"
    publishDir "${params.outdir}/enrichment", mode: 'copy'

    input:
    path gene_lists
    path universe
    path gene_sets
    val fdr
    val min_size
    val max_size
    val min_overlap

    output:
    path "custom_enrichment.tsv", emit: results

    script:
    """
    enrich_custom.py \\
        --gene_lists ${gene_lists} \\
        --universe ${universe} \\
        --gene_sets ${gene_sets} \\
        --fdr ${fdr} \\
        --min_size ${min_size} \\
        --max_size ${max_size} \\
        --min_overlap ${min_overlap}
    """
}

process ENRICH_HUMAN {
    tag "GO/KEGG"
    publishDir "${params.outdir}/enrichment", mode: 'copy'

    input:
    path gene_lists
    path universe
    val gene_id_type
    val run_kegg
    val fdr
    val min_size
    val max_size
    val min_overlap

    output:
    path "GO.tsv", emit: go
    path "KEGG.tsv", emit: kegg

    script:
    """
    enrich_human.R \\
        --gene_lists ${gene_lists} \\
        --universe ${universe} \\
        --gene_id_type ${gene_id_type} \\
        --run_kegg ${run_kegg} \\
        --fdr ${fdr} \\
        --min_size ${min_size} \\
        --max_size ${max_size} \\
        --min_overlap ${min_overlap}
    """
}
