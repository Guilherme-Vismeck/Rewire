process PLOT_NETWORK {
    tag "network_plot"
    publishDir "${params.outdir}/figures", mode: 'copy'

    input:
    path edges
    path stable
    path de_table
    val label_a
    val label_b
    val max_genes

    output:
    path "rewiring_network.png", emit: png
    path "rewiring_network.svg", emit: svg

    script:
    """
    plot_network.py \\
        --edges ${edges} \\
        --stable ${stable} \\
        --de_vs_rewiring ${de_table} \\
        --label_a ${label_a} \\
        --label_b ${label_b} \\
        --max_genes ${max_genes}
    """
}
