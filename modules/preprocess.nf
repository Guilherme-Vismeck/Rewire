process PREPROCESS {
    tag "preprocess"
    publishDir "${params.outdir}/preprocessing", mode: 'copy', pattern: 'filtered_expression.tsv'

    input:
    path expression
    path metadata
    path qc_report
    val group_col
    val group_a
    val group_b
    val top_variable_genes

    output:
    path "filtered_expression.tsv", emit: filtered
    path "${group_a}_expression.tsv", emit: expr_a
    path "${group_b}_expression.tsv", emit: expr_b

    script:
    """
    preprocess.py \\
        --expression ${expression} \\
        --metadata ${metadata} \\
        --group_col ${group_col} \\
        --group_a ${group_a} \\
        --group_b ${group_b} \\
        --top_variable_genes ${top_variable_genes}
    """
}
