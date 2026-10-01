process BUILD_NETWORK {
    tag "${label}"
    publishDir "${params.outdir}/networks", mode: 'copy', pattern: '*_network.tsv'

    input:
    tuple val(label), path(expression)
    val method
    val min_cor
    val fdr

    output:
    tuple val(label), path("${label}_correlations.tsv.gz"), emit: correlations
    path "${label}_network.tsv", emit: network

    script:
    """
    build_network.py \\
        --expression ${expression} \\
        --label ${label} \\
        --method ${method} \\
        --min_cor ${min_cor} \\
        --fdr ${fdr}
    """
}
