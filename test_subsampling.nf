#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

import groovy.json.JsonSlurper

include { SEQTK_SAMPLE } from './modules/local/seqtk_sample/main'
include { CHECK_PAIRING } from './modules/local/check_pairing/main'
include { COVERAGE_ANALYSIS } from './modules/local/coverage_analysis/main'
include { SRATOOLS_PREFETCH } from './modules/nf-core/sratools/prefetch/main'
include { SRATOOLS_FASTERQDUMP } from './modules/nf-core/sratools/fasterqdump/main'
include { FASTP } from './modules/nf-core/fastp/main'
include { BWA_MEM } from './modules/nf-core/bwa/mem/main'
include { IVAR_VARIANTS } from './modules/nf-core/ivar/variants/main'

def getFastaHeader(path) {
    return file(path).readLines().find { it.startsWith('>') }.substring(1).split()[0]
}

def getSraAccessions(jsonlPath) {
    def jsonSlurper = new JsonSlurper()
    def sraAccessions = []
    
    file(jsonlPath).readLines().each { line ->
        def json = jsonSlurper.parseText(line)
        if (json.sraAccessions) {
            sraAccessions.addAll(json.sraAccessions)
        }
    }
    
    return sraAccessions
}

workflow {
    def test_sras = getSraAccessions(params.sra_metadata)
    def subsample_levels = params.subsample_levels
    
    reference_ch = Channel.fromPath(params.reference)
    reference_val = Channel.value(params.reference)
    gff_val = Channel.value(params.ivar_gff)
    reference_header_val = Channel.value(getFastaHeader(params.reference))
    
    // download and process the test SRAs
    sra_input = Channel.from(test_sras)
        .map { sra -> [[id: sra, sra_accession: sra], sra] }
    
    SRATOOLS_PREFETCH(sra_input)
    SRATOOLS_FASTERQDUMP(SRATOOLS_PREFETCH.out.sra)
    
    FASTP(
        SRATOOLS_FASTERQDUMP.out.reads
            .map { meta, reads ->
                tuple(
                    meta,
                    reads,
                    params.fastp_trim_front_read_1,
                    params.fastp_trim_front_read_2,
                    params.fastp_trim_tail_read_1,
                    params.fastp_trim_tail_read_2
                )
            }
    )
    
    // Check pairing and fix improperly paired files
    CHECK_PAIRING(FASTP.out.trimmed_reads)
    
    // create subsampling combos
    trimmed_reads = CHECK_PAIRING.out.reads
    
    subsample_input = trimmed_reads
        .combine(Channel.from(subsample_levels))
        .map { meta, reads, num_reads ->
            def new_meta = meta.clone()
            new_meta.subsample = num_reads
            new_meta.id = "${meta.id}_${num_reads}"
            [new_meta, reads, num_reads]
        }
    
    // subsample reads
    subsample_split = subsample_input
        .multiMap { meta, reads, num_reads ->
            reads: [meta, reads]
            num_reads: num_reads
        }
    
    SEQTK_SAMPLE(
        subsample_split.reads,
        subsample_split.num_reads
    )
    
    // Run BWA alignment on subsampled reads
    BWA_MEM(
        SEQTK_SAMPLE.out.reads,
        reference_val
    )
    
    // variant calling prep
    bam_for_variants = BWA_MEM.out.bam
        .combine(reference_header_val)
        .map { meta, bam_file, header ->
            def gene = "genome"
            [meta, bam_file, gene, header]
        }
        .combine(gff_val)
        .map { meta, bam_file, gene, header, gff_file ->
            [meta, bam_file, gene, header, gff_file]
        }
    
    // Call variants w/ivar
    IVAR_VARIANTS(
        bam_for_variants,
        reference_val,
        params.variant_threshold,
        params.variant_min_depth
    )
    
    // Collect all variant files and run coverage analysis
    all_variants = IVAR_VARIANTS.out.variants.collect()
    
    COVERAGE_ANALYSIS(all_variants)
    
    IVAR_VARIANTS.out.variants.view()
} 