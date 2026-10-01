process DIFFERENTIAL_NETWORK {
    tag "differential"
    publishDir "${params.outdir}/differential_network", mode: 'copy'

    input:
    path corr_a
    path corr_b
    path expr_a
    path expr_b
    val label_a
    val label_b
    val fdr

    output:
    path "differential_edges.tsv", emit: edges
    path "gained_edges.tsv"
    path "lost_edges.tsv"
    path "sign_flips.tsv"
    path "gene_rewiring.tsv", emit: genes

    script:
    """
    N_A=\$(head -n 1 ${expr_a} | awk '{print NF-1}')
    N_B=\$(head -n 1 ${expr_b} | awk '{print NF-1}')

    differential_correlation.py \\
        --corr_a ${corr_a} \\
        --corr_b ${corr_b} \\
        --label_a ${label_a} \\
        --label_b ${label_b} \\
        --n_a \$N_A \\
        --n_b \$N_B \\
        --fdr ${fdr}
    """
}
