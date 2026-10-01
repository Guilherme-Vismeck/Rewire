process VALIDATE_INPUTS {
    tag "validate"
    publishDir "${params.outdir}/qc", mode: 'copy'

    input:
    path expression
    path metadata
    val group_col
    val group_a
    val group_b

    output:
    path "qc_report.tsv", emit: report

    script:
    """
    validate_inputs.py \\
        --expression ${expression} \\
        --metadata ${metadata} \\
        --group_col ${group_col} \\
        --group_a ${group_a} \\
        --group_b ${group_b} \\
        --report qc_report.tsv
    """
}
