include { VALIDATE_INPUTS } from '../modules/validate_inputs.nf'
include { PREPROCESS } from '../modules/preprocess.nf'
include { BUILD_NETWORK } from '../modules/build_network.nf'
include { DIFFERENTIAL_NETWORK } from '../modules/differential_correlation.nf'
include { TOPOLOGY } from '../modules/topology.nf'
include { COMMUNITIES } from '../modules/communities.nf'
include { BOOTSTRAP } from '../modules/bootstrap.nf'
include { DIFFERENTIAL_EXPRESSION } from '../modules/differential_expression.nf'
include { INTEGRATE_DE_RW } from '../modules/integrate_de_rewiring.nf'
include { PREPARE_GENE_LISTS; ENRICH_CUSTOM; ENRICH_HUMAN } from '../modules/enrichment.nf'
include { EXPORT_NETWORK } from '../modules/export_network.nf'
include { REPORT } from '../modules/report.nf'

workflow REWIRENF {
    if( !params.expression || !params.metadata || !params.group_a || !params.group_b ) {
        error "Parâmetros obrigatórios: --expression, --metadata, --group_a, --group_b"
    }
    if( params.run_enrichment ) {
        if( !params.run_de ) {
            error "--run_enrichment requer --run_de true (as listas usam as categorias DE x rewiring)"
        }
        if( params.organism && params.organism != 'human' ) {
            error "Organismo '${params.organism}' sem banco embutido. Use --gene_sets com seus próprios conjuntos (colunas term e gene, ou .gmt)"
        }
        if( !params.gene_sets && params.organism != 'human' ) {
            error "--run_enrichment requer --organism human e/ou --gene_sets <arquivo>"
        }
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
        .filter { label, _f -> label == params.group_a }
        .map { _label, f -> f }
    def corr_b = BUILD_NETWORK.out.correlations
        .filter { label, _f -> label == params.group_b }
        .map { _label, f -> f }

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

        if( params.run_enrichment ) {
            PREPARE_GENE_LISTS(
                INTEGRATE_DE_RW.out.table,
                BOOTSTRAP.out.stable,
                COMMUNITIES.out.gene_modules
            )

            if( params.gene_sets ) {
                ENRICH_CUSTOM(
                    PREPARE_GENE_LISTS.out.gene_lists,
                    PREPARE_GENE_LISTS.out.universe,
                    file(params.gene_sets, checkIfExists: true),
                    params.enrich_fdr,
                    params.enrich_min_size,
                    params.enrich_max_size,
                    params.enrich_min_overlap
                )
            }

            if( params.organism == 'human' ) {
                ENRICH_HUMAN(
                    PREPARE_GENE_LISTS.out.gene_lists,
                    PREPARE_GENE_LISTS.out.universe,
                    params.gene_id_type,
                    params.run_kegg,
                    params.enrich_fdr,
                    params.enrich_min_size,
                    params.enrich_max_size,
                    params.enrich_min_overlap
                )
            }
        }
    }

    EXPORT_NETWORK(
        net_a,
        net_b,
        DIFFERENTIAL_NETWORK.out.edges,
        DIFFERENTIAL_NETWORK.out.genes,
        params.group_a,
        params.group_b
    )

    if( params.run_report ) {
        def optional_de = params.run_de
            ? INTEGRATE_DE_RW.out.table.mix(INTEGRATE_DE_RW.out.summary, DIFFERENTIAL_EXPRESSION.out.results)
            : channel.empty()
        def optional_custom = (params.run_de && params.run_enrichment && params.gene_sets)
            ? ENRICH_CUSTOM.out.results
            : channel.empty()
        def optional_human = (params.run_de && params.run_enrichment && params.organism == 'human')
            ? ENRICH_HUMAN.out.go.mix(ENRICH_HUMAN.out.kegg)
            : channel.empty()

        def report_inputs = VALIDATE_INPUTS.out.report
            .mix(
                DIFFERENTIAL_NETWORK.out.edges,
                DIFFERENTIAL_NETWORK.out.genes,
                TOPOLOGY.out.metrics,
                TOPOLOGY.out.turnover,
                BOOTSTRAP.out.edge_stability,
                BOOTSTRAP.out.gene_stability,
                BOOTSTRAP.out.stable,
                COMMUNITIES.out.gene_modules,
                COMMUNITIES.out.module_summary,
                BUILD_NETWORK.out.network,
                BUILD_NETWORK.out.correlations.map { _label, f -> f },
                optional_de,
                optional_custom,
                optional_human
            )
            .collect()

        def params_text = params.findAll { _k, v -> v != null }
            .collect { k, v -> "${k}\t${v}" }
            .sort()
            .join('\n')
        def params_file = channel.of(params_text).collectFile(name: 'parameters.tsv', newLine: true)

        REPORT(
            report_inputs,
            params_file,
            params.group_a,
            params.group_b,
            workflow.manifest.version,
            workflow.nextflow.version
        )
    }
}
