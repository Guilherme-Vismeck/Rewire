process DIFFERENTIAL_EXPRESSION {
    tag "limma"
    publishDir "${params.outdir}/differential_expression", mode: 'copy'

    input:
    path expression
    path metadata
    val group_col
    val group_a
    val group_b
    val de_fdr
    val min_logfc

    output:
    path "DE_results.tsv", emit: results

    script:
    """
    differential_expression.R \\
        --expression ${expression} \\
        --metadata ${metadata} \\
        --group_col ${group_col} \\
        --group_a ${group_a} \\
        --group_b ${group_b} \\
        --de_fdr ${de_fdr} \\
        --min_logfc ${min_logfc}
    """
}
