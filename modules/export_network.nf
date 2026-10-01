process EXPORT_NETWORK {
    tag "graphml"
    publishDir "${params.outdir}/networks_graphml", mode: 'copy'

    input:
    path network_a
    path network_b
    path edges
    path genes
    val label_a
    val label_b

    output:
    path "*.graphml"

    script:
    """
    export_network.py \\
        --network_a ${network_a} \\
        --network_b ${network_b} \\
        --edges ${edges} \\
        --genes ${genes} \\
        --label_a ${label_a} \\
        --label_b ${label_b}
    """
}
