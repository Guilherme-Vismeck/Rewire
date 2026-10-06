process REPORT {
    tag "report"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path inputs, stageAs: 'inputs/*'
    path parameters
    val label_a
    val label_b
    val pipeline_version
    val nextflow_version

    output:
    path "report.html", emit: report

    script:
    """
    build_report.py \\
        --inputs inputs \\
        --parameters ${parameters} \\
        --label_a ${label_a} \\
        --label_b ${label_b} \\
        --pipeline_version ${pipeline_version} \\
        --nextflow_version ${nextflow_version} \\
        --output report.html
    """
}
