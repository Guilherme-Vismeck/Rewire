process TOPOLOGY {
    tag "topology"
    publishDir "${params.outdir}/topology", mode: 'copy', pattern: 'network_metrics.tsv'
    publishDir "${params.outdir}/rewiring", mode: 'copy', pattern: 'neighbor_turnover.tsv'

    input:
    path network_a
    path network_b
    path expression
    val label_a
    val label_b

    output:
    path "network_metrics.tsv", emit: metrics
    path "neighbor_turnover.tsv", emit: turnover

    script:
    """
    topology.py \\
        --network_a ${network_a} \\
        --network_b ${network_b} \\
        --expression ${expression} \\
        --label_a ${label_a} \\
        --label_b ${label_b}
    """
}
