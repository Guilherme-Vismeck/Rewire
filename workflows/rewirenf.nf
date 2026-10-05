include { VALIDATE_INPUTS } from '../modules/validate_inputs.nf'
include { PREPROCESS } from '../modules/preprocess.nf'
include { BUILD_NETWORK } from '../modules/build_network.nf'
include { DIFFERENTIAL_NETWORK } from '../modules/differential_correlation.nf'
include { TOPOLOGY } from '../modules/topology.nf'
include { COMMUNITIES } from '../modules/communities.nf'
include { BOOTSTRAP } from '../modules/bootstrap.nf'
include { DIFFERENTIAL_EXPRESSION } from '../modules/differential_expression.nf'
include { INTEGRATE_DE_RW } from '../modules/integrate_de_rewiring.nf'
include { EXPORT_NETWORK } from '../modules/export_network.nf'

workflow REWIRENF {
    if( !params.expression || !params.metadata || !params.group_a || !params.group_b ) {
        error "Parâmetros obrigatórios: --expression, --metadata, --group_a, --group_b"
    }

    def expression = file(params.expression, checkIfExists: true)
    def metadata = file(params.metadata, checkIfExists: true)

    VALIDATE_INPUTS(expression, metadata, params.group_col, params.group_a, params.group_b)

    PREPROCESS(
        expression,
        metadata,
        VALIDATE_INPUTS.out.report,
        params.group_col,
        params.group_a,
        params.group_b,
        params.top_variable_genes
    )

    def groups = PREPROCESS.out.expr_a
        .map { f -> tuple(params.group_a, f) }
        .mix( PREPROCESS.out.expr_b.map { f -> tuple(params.group_b, f) } )

    BUILD_NETWORK(groups, params.method, params.min_cor, params.fdr)

    def corr_a = BUILD_NETWORK.out.correlations
        .filter { label, f -> label == params.group_a }
        .map { label, f -> f }
    def corr_b = BUILD_NETWORK.out.correlations
        .filter { label, f -> label == params.group_b }
        .map { label, f -> f }

    DIFFERENTIAL_NETWORK(
        corr_a,
        corr_b,
        PREPROCESS.out.expr_a,
        PREPROCESS.out.expr_b,
        params.group_a,
        params.group_b,
        params.fdr
    )

    def net_a = BUILD_NETWORK.out.network.filter { f -> f.name == "${params.group_a}_network.tsv" }
    def net_b = BUILD_NETWORK.out.network.filter { f -> f.name == "${params.group_b}_network.tsv" }

    TOPOLOGY(
        net_a,
        net_b,
        PREPROCESS.out.filtered,
        params.group_a,
        params.group_b
    )

    BOOTSTRAP(
        DIFFERENTIAL_NETWORK.out.edges,
        PREPROCESS.out.expr_a,
        PREPROCESS.out.expr_b,
        params.group_a,
        params.group_b,
        params.method,
        params.min_cor,
        params.bootstrap,
        params.bootstrap_seed,
        params.stability_threshold
    )

    if( params.run_de ) {
        DIFFERENTIAL_EXPRESSION(
            PREPROCESS.out.filtered,
            metadata,
            params.group_col,
            params.group_a,
            params.group_b,
            params.de_fdr,
            params.min_logfc
        )

        INTEGRATE_DE_RW(
            DIFFERENTIAL_EXPRESSION.out.results,
            BOOTSTRAP.out.gene_stability,
            params.min_rewired_edges
        )
    }

    COMMUNITIES(
        net_a,
        net_b,
        PREPROCESS.out.filtered,
        params.group_a,
        params.group_b,
        params.community_resolution,
        params.community_seed,
        params.module_overlap_threshold
    )

    EXPORT_NETWORK(
        net_a,
        net_b,
        DIFFERENTIAL_NETWORK.out.edges,
        DIFFERENTIAL_NETWORK.out.genes,
        params.group_a,
        params.group_b
    )
}
