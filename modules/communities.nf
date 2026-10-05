process COMMUNITIES {
    tag "communities"
    publishDir "${params.outdir}/communities", mode: 'copy'

    input:
    path network_a
    path network_b
    path expression
    val label_a
    val label_b
    val resolution
    val seed
    val overlap_threshold

    output:
    path "gene_modules.tsv", emit: gene_modules
    path "module_summary.tsv", emit: module_summary

    script:
    """
    communities.py \\
        --network_a ${network_a} \\
        --network_b ${network_b} \\
        --expression ${expression} \\
        --label_a ${label_a} \\
        --label_b ${label_b} \\
        --resolution ${resolution} \\
        --seed ${seed} \\
        --module_overlap_threshold ${overlap_threshold}
    """
}
