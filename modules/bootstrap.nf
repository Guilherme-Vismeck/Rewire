process BOOTSTRAP {
    tag "bootstrap"
    publishDir "${params.outdir}/bootstrap", mode: 'copy', pattern: '*_stability.tsv'
    publishDir "${params.outdir}/rewiring", mode: 'copy', pattern: 'stable_rewiring.tsv'

    input:
    path edges
    path expr_a
    path expr_b
    val label_a
    val label_b
    val method
    val min_cor
    val n_boot
    val seed
    val threshold

    output:
    path "edge_stability.tsv", emit: edge_stability
    path "gene_stability.tsv", emit: gene_stability
    path "stable_rewiring.tsv", emit: stable

    script:
    """
    bootstrap.py \\
        --edges ${edges} \\
        --expr_a ${expr_a} \\
        --expr_b ${expr_b} \\
        --label_a ${label_a} \\
        --label_b ${label_b} \\
        --method ${method} \\
        --min_cor ${min_cor} \\
        --bootstrap ${n_boot} \\
        --seed ${seed} \\
        --stability_threshold ${threshold}
    """
}
