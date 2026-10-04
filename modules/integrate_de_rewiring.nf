process INTEGRATE_DE_RW {
    tag "integrate"
    publishDir "${params.outdir}/rewiring", mode: 'copy'

    input:
    path de_results
    path gene_stability
    val min_rewired_edges

    output:
    path "de_vs_rewiring.tsv", emit: table
    path "category_summary.tsv", emit: summary

    script:
    """
    integrate_de_rewiring.py \\
        --de ${de_results} \\
        --gene_stability ${gene_stability} \\
        --min_rewired_edges ${min_rewired_edges}
    """
}
