/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PLASMID ANALYSIS SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MOB-suite reconstruction -> per-sample contig/mobtyper reports -> summary
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { MOB_RECON       } from '../../../modules/local/mob_recon'
include { PLASMID_SUMMARY } from '../../../modules/local/plasmid_summary'

workflow PLASMID_ANALYSIS {

    take:
    ch_assemblies  // channel: [ val(meta), path(assembly) ]

    main:

    MOB_RECON ( ch_assemblies )

    PLASMID_SUMMARY (
        MOB_RECON.out.full_output.map { meta, dir -> dir }.collect()
    )

    emit:
    mob_output      = MOB_RECON.out.full_output     // [ val(meta), path(dir) ]
    contig_report   = MOB_RECON.out.contig_report   // [ val(meta), path(txt) ]
    mobtyper        = MOB_RECON.out.mobtyper         // [ val(meta), path(txt) ]
    plasmid_summary = PLASMID_SUMMARY.out.summary   // path(plasmid_summary.tsv)
}
