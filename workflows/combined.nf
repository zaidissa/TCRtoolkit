/*
 * COMBINED workflow (Scenario 3)
 *
 * Both bulk and single-cell TCR data are available and run in parallel.
 *
 * Track A (Bulk):  INPUT_CHECK -> CONVERT -> ANNOTATE -> SAMPLE -> PATIENT -> COMPARE
 * Track B (SC):    VDJ_QC -> TCELL_INTEGRATION
 *                    -> Bridge 1 (SC_TO_BULK) -> ANNOTATE_SC -> PATIENT_SC (GIANA+GLIPH2)
 *                                            -> TCRDIST3_MATRIX
 *                    -> Bridge 2+ (CLUSTER_TO_SC) -> enriched Seurat
 *                    -> CONGA -> CONSENSUS -> REPERTOIRE
 *
 * Merge point: MASTER_SUMMARY collects outputs from both tracks.
 *
 * Inputs:
 *   Bulk  - --samplesheet, --input_format
 *   SC    - --input_annotated_object, --input_vdj_contigs, --sample_sheet
 */

include { INPUT_CHECK }          from '../subworkflows/bulk/input_check.nf'
include { CONVERT }              from '../subworkflows/bulk/convert.nf'
include { ANNOTATE as ANNOTATE_BULK } from '../subworkflows/bulk/annotate.nf'
include { ANNOTATE as ANNOTATE_SC }   from '../subworkflows/bulk/annotate.nf'
include { SAMPLE }               from '../subworkflows/bulk/sample.nf'
include { PATIENT as PATIENT_BULK }   from '../subworkflows/bulk/patient.nf'
include { PATIENT as PATIENT_SC }     from '../subworkflows/bulk/patient.nf'
include { COMPARE }              from '../subworkflows/bulk/compare.nf'
include { PSEUDOBULK_PHENOTYPE } from '../subworkflows/bulk/pseudobulk_phenotype.nf'
include { PSEUDOBULK_QC_SW }     from '../subworkflows/bulk/pseudobulk_qc.nf'

include { VDJ_QC_SW }            from '../subworkflows/scratch/vdj_qc.nf'
include { TCELL_INTEGRATION_SW } from '../subworkflows/scratch/tcell_integration.nf'
include { CONGA_SW }             from '../subworkflows/scratch/conga.nf'
include { CONSENSUS_SW }         from '../subworkflows/scratch/consensus_clustering.nf'
include { REPERTOIRE_SW }        from '../subworkflows/scratch/repertoire.nf'
include { MASTER_SUMMARY_SW }    from '../subworkflows/scratch/master_summary.nf'

include { SC_TO_BULK_SW }        from '../subworkflows/bridges/sc_to_bulk.nf'
include { CLUSTER_TO_SC_SW }     from '../subworkflows/bridges/cluster_to_sc.nf'

include { TCRDIST3_MATRIX }      from '../modules/bulk/sample/tcrdist3'

workflow COMBINED_WORKFLOW {

    def enabled     = { x -> x == null || x == true }
    def nofile      = file("${projectDir}/assets/NO_FILE")
    def input_format = params.input_format.toLowerCase()

    // ── Validate mandatory inputs ─────────────────────────────────────────
    if (!params.samplesheet)             error "Combined mode: provide --samplesheet (bulk)"
    if (!params.input_annotated_object)  error "Combined mode: provide --input_annotated_object (SC)"
    if (!params.input_vdj_contigs)       error "Combined mode: provide --input_vdj_contigs (SC)"
    if (!params.sample_sheet)            error "Combined mode: provide --sample_sheet (SC)"

    // ══════════════════════════════════════════════════════════════════════
    // TRACK A - Bulk TCRtoolkit analysis
    // ══════════════════════════════════════════════════════════════════════

    def levels = params.workflow_level.toLowerCase().tokenize(',')

    INPUT_CHECK( file(params.samplesheet) )

    if (input_format == 'adaptive' || input_format == 'cellranger') {
        CONVERT( INPUT_CHECK.out.sample_map, input_format )
        bulk_sample_map = CONVERT.out.sample_map_converted

        if (input_format == 'cellranger' && params.sobject_gex) {
            PSEUDOBULK_PHENOTYPE(
                CONVERT.out.pseudobulk_phenotype_files,
                INPUT_CHECK.out.samplesheet_utf8,
                levels
            )
        }
    } else {
        bulk_sample_map = INPUT_CHECK.out.sample_map
    }

    ANNOTATE_BULK( bulk_sample_map )

    if (levels.contains('sample')) {
        SAMPLE( bulk_sample_map, ANNOTATE_BULK.out.cdr3_pgen, ANNOTATE_BULK.out.olga_stats )
    }

    if (levels.contains('patient')) {
        PATIENT_BULK( ANNOTATE_BULK.out.processed_samples )
    }

    if (levels.contains('compare')) {
        COMPARE( ANNOTATE_BULK.out.concat_cdr3_sorted, ANNOTATE_BULK.out.cdr3_pgen )
    }

    // ══════════════════════════════════════════════════════════════════════
    // TRACK B - SCRATCH single-cell analysis (parallel with Track A)
    // ══════════════════════════════════════════════════════════════════════

    ch_annotated_object = Channel.fromPath(params.input_annotated_object, checkIfExists: true)
    ch_sample_sheet     = Channel.fromPath(params.sample_sheet,           checkIfExists: true)
    ch_project_name     = Channel.value(params.project_name)

    // Step B-1: VDJ QC
    vdj_qc_out = VDJ_QC_SW(
        ch_sample_sheet,
        ch_project_name,
        ch_annotated_object
    )

    def vdj_qc_per_sample_compact         = vdj_qc_out.qc_tables.flatten().filter { it.name == 'vdj_qc_per_sample_compact.tsv'      }.ifEmpty(nofile)
    def vdj_qc_before_after_summary       = vdj_qc_out.qc_tables.flatten().filter { it.name == 'qc_contigs_before_after_summary.tsv'}.ifEmpty(nofile)
    def vdj_qc_sample_sheet_resolved      = vdj_qc_out.qc_tables.flatten().filter { it.name == 'sample_sheet_resolved.tsv'          }.ifEmpty(nofile)
    def vdj_qc_clone_rank_abundance       = vdj_qc_out.qc_tables.flatten().filter { it.name == 'clone_rank_abundance.tsv'           }.ifEmpty(nofile)
    def vdj_qc_before_after_retention_fig = vdj_qc_out.qc_figures.flatten().filter { it.name == 'qc_before_after_retention.png'     }.ifEmpty(nofile)
    def vdj_qc_pairing_bar_fig            = vdj_qc_out.qc_figures.flatten().filter { it.name == 'pairing_bar_by_sample.png'         }.ifEmpty(nofile)
    def vdj_qc_clone_rank_abundance_fig   = vdj_qc_out.qc_figures.flatten().filter { it.name == 'clone_rank_abundance.png'          }.ifEmpty(nofile)
    def vdj_qc_multiple_chains_fig        = vdj_qc_out.qc_figures.flatten().filter { it.name == 'multiple_chains_by_sample.png'     }.ifEmpty(nofile)

    // Step B-2: T-cell integration
    tcell_out = TCELL_INTEGRATION_SW(
        vdj_qc_out.contigs_after_qc,
        ch_annotated_object,
        ch_project_name
    )

    // Step B-3: Bridge 1 - SC export_cells -> per-sample AIRR TSVs
    SC_TO_BULK_SW( tcell_out.export_cells )

    // Step B-3b: TCRtoolkit QC gate on pseudobulk-derived bulk TCR
    PSEUDOBULK_QC_SW( SC_TO_BULK_SW.out.sample_map )

    // Step B-4: TCRtoolkit bulk clustering on QC-passed SC-derived AIRR files
    ANNOTATE_SC( PSEUDOBULK_QC_SW.out.sample_map )

    PATIENT_SC( ANNOTATE_SC.out.processed_samples )

    TCRDIST3_MATRIX(
        ANNOTATE_SC.out.processed_samples,
        params.matrix_sparsity,
        params.distance_metric,
        file(params.db_path)
    )

    // Step B-5: Bridge 2+ - map cluster results -> enriched Seurat
    CLUSTER_TO_SC_SW(
        tcell_out.seurat_tcells_with_tcr,
        tcell_out.export_cells,
        PATIENT_SC.out.giana_clusters,
        PATIENT_SC.out.gliph2_cluster_details,
        TCRDIST3_MATRIX.out.clone_df,
        TCRDIST3_MATRIX.out.tcrdist_output.map { _meta, f -> f }
    )

    enriched_seurat = CLUSTER_TO_SC_SW.out.enriched_seurat

    // Step B-6: CoNGA on enriched Seurat
    conga_done = Channel.empty()
    if (enabled(params.run_conga)) {
        conga_out  = CONGA_SW( enriched_seurat, tcell_out.export_cells, ch_project_name )
        conga_done = conga_out.report_html
    }

    // Step B-7: Consensus clustering
    consensus_done = Channel.empty()
    if (enabled(params.run_consensus)) {
        def gliph2_export  = CLUSTER_TO_SC_SW.out.gliph2_export.ifEmpty(nofile)
        def tcrdist_export = CLUSTER_TO_SC_SW.out.tcrdist_export.ifEmpty(nofile)
        def giana_export   = CLUSTER_TO_SC_SW.out.giana_export.ifEmpty(nofile)

        CONSENSUS_SW(
            enriched_seurat,
            tcell_out.export_cells,
            gliph2_export,
            tcrdist_export,
            giana_export,
            ch_project_name
        )
        consensus_done = CONSENSUS_SW.out.report_html
    }

    downstream_seurat = enabled(params.run_consensus)
        ? CONSENSUS_SW.out.seurat_with_consensus
        : enriched_seurat

    downstream_export = enabled(params.run_consensus)
        ? CONSENSUS_SW.out.export_cells
        : tcell_out.export_cells

    // Step B-8: Repertoire
    repertoire_done = Channel.empty()
    if (enabled(params.run_repertoire)) {
        REPERTOIRE_SW( downstream_seurat, downstream_export, ch_project_name )
        repertoire_done = REPERTOIRE_SW.out.report_html
    }

    // ══════════════════════════════════════════════════════════════════════
    // MERGE - Master summary aggregates both tracks
    // ══════════════════════════════════════════════════════════════════════
    if (enabled(params.run_master_summary)) {
        master_barrier = conga_done
            .mix(consensus_done)
            .mix(repertoire_done)
            .collect()

        MASTER_SUMMARY_SW(
            downstream_seurat,
            downstream_export,
            vdj_qc_per_sample_compact,
            vdj_qc_before_after_summary,
            vdj_qc_sample_sheet_resolved,
            vdj_qc_clone_rank_abundance,
            vdj_qc_before_after_retention_fig,
            vdj_qc_pairing_bar_fig,
            vdj_qc_clone_rank_abundance_fig,
            vdj_qc_multiple_chains_fig,
            master_barrier,
            ch_project_name
        )
    }
}
