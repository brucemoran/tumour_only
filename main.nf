#!/usr/bin/env nextflow

def helpMessage() {
  log.info"""
  -----------------------------------------------------------------------
                          TUMOUR_ONLY PIPELINE
  -----------------------------------------------------------------------
  Usage:

  nextflow run brucemoran/tumour_only

  Mandatory arguments:

    -profile        [str]       Configuration profile
                                (required: standard,singularity)

    --sampleCsv     [file]      CSV format, headers: sampleID, meta, read1 (e.g.
                                /path/to/read1.fastq.gz), read2
                                (e.g. /path/to/read2.fastq.gz); use meta for
                                naming in PCGR, CPSR reports

    --runID         [str]       Name for run, used to tag outputs

    --refDir        [file]      Path of dir in which reference data are held;
                                this should be created by download-references.nf
                                and contain dir <assembly>

    --assembly      [str]       Either GRCh37 or GRCh38 (default), as per
                                download-references.nf

    --outDir        [str]       Path to send outputs (default: ./)

    --email         [str]       Email address to send reports

  General Optional Arguments:

    --scatGath      [int]       Number of pieces to divide intervalList into for
                                scattering to variant calling processes
                                (default: 20 for exome/panel, 100 for WGS)

    --sampleCat     [str]       File used when fastq data is in multiple files
                                which are cat'ed; replaces --sampleCsv; headers:
                                sampleID, meta, dir (contains fastq to be cat'ed),
                                ext (extension scheme for parsing read1, read2 e.g.
                                _1.fq.gz;_2.fq.gz).

    --multiqcConfig [str]       Config file for multiqc
                                (default: bin/tumour_only.multiQC_config.yaml)

    --seqLevel      [str]       WGS (default), exome or panel; used to tag output
                                and for specifying which reference data to use

    --levelTag      [str]       Tag used for exome/panel kit when running
                                download-references.nf

    """.stripIndent()
}

if (params.help) exit 0, helpMessage()

//Test Mandatory Arguments
if(params.sampleCsv && params.sampleCat){
  exit 1, "Please include only one of --sampleCsv or --sampleCat, see --help for format"
}

if(params.sampleCsv == null && params.sampleCat == null){
  exit 1, "Please include one of --sampleCsv or --sampleCat, see --help for format"
}

if(!Channel.from(params.runID, checkIfExists: true)){
    exit 1, "Please include --runID <your_runID>"
}

if(!Channel.from(params.refDir, checkIfExists: true)){
  exit 1, "Please include --refDir <path> see github.com/brucemoran/tumour_only/ for how to run download-references.nf"
}

if(!Channel.from(params.assembly, checkIfExists: true)){
    exit 1, "Please include --assembly <GRCh3x>"
}

if(!params.email){
    exit 1, "Please include --email your@email.com"
}

if(params.seqLevel != "WGS" && params.levelTag == null){
    exit 1, "Please define --levelTag when using --seqLevel exome/panel"
}
//Global Variables based on input
params.outdir = "${params.outDir}/${params.seqLevel}_output"
params.seqlevel = "${params.seqLevel}".toLowerCase()

//Java task memory allocation via task.memory
javaTaskmem = { it.replace(" GB", "g") }

//Reference data as value channels and reusable therefore
reference = [
    grchvers: false,
    fa: false,
    fai: false,
    dict: false,
    bwa: false,
    hc_dbs: false,
    dbsnp: false,
    pcgrbase: false,
    intlist: false,
    seqlevel: false,
    bbres: false
]

reference.grchvers  = Channel.fromPath("${params.refDir}/${params.assembly}/pcgr/data/*", type: 'dir').getVal()
reference.fa = Channel.value(file(params.genomes[params.assembly].fa))
reference.fai = Channel.value(file(params.genomes[params.assembly].fai))
reference.dict = Channel.value(file(params.genomes[params.assembly].dict))
reference.bwa = Channel.value(file(params.genomes[params.assembly].bwa))
reference.hc_dbs = Channel.value(file(params.genomes[params.assembly].hc_dbs))
reference.dbsnp = Channel.value(file(params.genomes[params.assembly].dbsnp))
reference.pcgrbase = Channel.value(file(params.genomes[params.assembly].pcgr))
reference.pathseq = Channel.value(file(params.genomes[params.assembly].pathseq))
reference.refflat = Channel.value(file(params.genomes[params.assembly].refflat))

//if seqlevel is not WGS, there is a params.genomes dir per exome or panel
//holding a levelTag dir which we need to specify
reference.seqlevel = Channel.value(file(params.genomes[params.assembly]."${params.seqlevel}"))

//set cosmic
reference.cosmic = params.cosmic == true ? Channel.value(file(params.genomes[params.assembly].cosmic)) : null

//setting of intlist, bed based on seqlevel and levelTag
reference.intlist = params.seqlevel == "wgs" ? Channel.fromPath("${params.refDir}/${params.assembly}/${params.seqlevel}/wgs.bed.interval_list").getVal() : Channel.fromPath("${params.refDir}/${params.assembly}/${params.seqlevel}/${params.levelTag}/${params.levelTag}.bed.interval_list").getVal()
reference.bed = params.seqlevel == "wgs" ? Channel.fromPath("${params.refDir}/${params.assembly}/${params.seqlevel}/wgs.bed").getVal() : Channel.fromPath("${params.refDir}/${params.assembly}/${params.seqlevel}/${params.levelTag}/${params.levelTag}.bed").getVal()

/*
================================================================================
                          -0. PREPROCESS INPUT SAMPLE FILE
================================================================================
*/
/* 0.00: Input using sample.csv, bam.csv, sample_cat
*/
if(params.sampleCsv){
  Channel.fromPath("${params.sampleCsv}")
         .splitCsv( header: true )
         .map { row -> [row.sampleID, row.meta, file(row.read1), file(row.read2)] }
         .set { bbduking }
}

if(params.sampleCat){
  Channel.fromPath("${params.sampleCat}")
         .splitCsv( header: true )
         .map { row -> [row.sampleID, row.meta, row.dir, row.ext] }
         .set { samplecating }

  process Samplecat {

    label 'low_mem'
    publishDir "${params.outdir}/samples/${sampleID}/cat", mode: "copy"

    input:
    tuple val(sampleID), val(meta), val(dir), val(ext) from samplecating

    output:
    tuple val(sampleID), val(meta), file(read1), file(read2) into bbduking

    script:
    rd1ext = "${ext}".split(';')[0]
    rd2ext = "${ext}".split(';')[1]
    read1 = "${sampleID}.R1.fastq.gz"
    read2 = "${sampleID}.R2.fastq.gz"
    """
    #! bash
    cat \$(find ${dir} | grep ${rd1ext} | sort) > ${read1}
    cat \$(find ${dir} | grep ${rd2ext} | sort) > ${read2}
    """
  }
}

/*
================================================================================
                          0. PREPROCESS INPUT FASTQ
================================================================================
*/
// 0.1: Input trimming
process bbduk {

  label 'med_mem'
  publishDir path: "${params.outdir}/samples/${sampleID}/bbduk", mode: "copy", pattern: "*.txt"

  input:
  tuple val(sampleID), val(meta), file(read1), file(read2) from bbduking

  output:
  file('*.txt') into log_bbduk
  tuple val(sampleID), val(meta), file('*.bbduk.R1.fastq.gz'), file('*.bbduk.R2.fastq.gz') into bwa_memming
  tuple val(sampleID), val(meta), file('*.bbduk.R1.fastq.gz'), file('*.bbduk.R2.fastq.gz'), file(read1), file(read2) into fastping
  tuple val(sampleID), val(meta), file(read1), file(read2) into fastqcing
  tuple val(sampleID), val(meta) into pcgr_meta

  script:
  def taskmem = task.memory == null ? "" : "-Xmx" + javaTaskmem("${task.memory}")
  """
  {
  sh bbduk.sh ${taskmem} \
    in1=${read1} \
    in2=${read2} \
    out1=${sampleID}".bbduk.R1.fastq.gz" \
    out2=${sampleID}".bbduk.R2.fastq.gz" \
    k=31 \
    mink=5 \
    hdist=1 \
    ktrim=r \
    trimq=20 \
    qtrim=rl \
    maq=20 \
    ref=/opt/miniconda/envs/somatic_n-of-1/opt/bbmap-adapters.fa \
    tpe \
    tbo \
    stats=${sampleID}".bbduk.adapterstats.txt" \
    overwrite=T
  } 2>&1 | tee > ${sampleID}.bbduk.runstats.txt
  """
}

// 0.2: fastp QC of pre-, post-bbduk
process fastp {

  label 'low_mem'
  publishDir "${params.outdir}/samples/${sampleID}/fastp", mode: "copy", pattern: "*.html"

  input:
  tuple val(sampleID), val(meta), file(preread1), file(preread2), file(postread1), file(postread2) from fastping

  output:
  file('*.html') into fastp_html
  file('*.json') into fastp_multiqc

  script:
  """
  fastp -w ${task.cpus} -h ${sampleID}"_pre.fastp.html" -j ${sampleID}"_pre.fastp.json" --in1 ${preread1} --in2 ${preread2}

  fastp -w ${task.cpus} -h ${sampleID}"_post.fastp.html" -j ${sampleID}"_post.fastp.json" --in1 ${postread1} --in2 ${postread2}
  """
}

// 0.3: fastQC of per, post-bbduk
process fastqc {

  label 'low_mem'
  publishDir "${params.outdir}/samples/${sampleID}/fastqc", mode: "copy", pattern: "*.html"

  input:
  tuple val(sampleID), val(meta), file(read1), file(read2) from fastqcing

  output:
  file('*.html') into fastqc_multiqc

  script:
  """
  #!/bin/bash
  fastqc -t ${task.cpus} --noextract -o ./ ${read1} ${read2}
  """
}

/*
================================================================================
                        1. ALIGNMENT AND BAM PROCESSING
================================================================================
*/
// 1.0: Input alignment
process bwamem {

  label 'high_mem'
  errorStrategy 'retry'
  maxRetries 3

  input:
  tuple val(sampleID), val(meta), file(read1), file(read2) from bwa_memming
  file(bwa) from reference.bwa

  output:
  tuple val(sampleID), val(meta), file('*.bam'), file('*.bai') into (cramming, dup_marking)

  script:
  def fa = "${bwa}/*fasta"
  """
  DATE=\$(date +"%Y-%m-%dT%T")
  RGLINE="@RG\\tID:${sampleID}\\tPL:ILLUMINA\\tSM:${sampleID}\\tDS:tumour_only\\tCN:UCD\\tLB:LANE_X\\tDT:\$DATE"

  bwa mem \
    -t${task.cpus} \
    -M \
    -R \$RGLINE \
    ${fa} \
    ${read1} ${read2} | \
    samtools sort -T "tmp."${sampleID} -o ${sampleID}".sort.bam"
  samtools index ${sampleID}".sort.bam"
  """
}

// 1.1: CRAM alignment and output
// TODO: output upload schema for ENA/EGA
process Cram {

  label 'low_mem'
  publishDir path: "${params.outdir}/samples/${sampleID}/bwa", mode: "copy", pattern: "*.cra*"

  input:
  tuple val(sampleID), val(meta), file(bam), file(bai) from cramming
  file(bwa) from reference.bwa

  output:
  tuple file('*.cram'), file('*.crai') into completedcram

  script:
  """
  samtools view -hC -T ${bwa}/*fasta ${sampleID}".sort.bam" > ${sampleID}".sort.cram"
  samtools index ${sampleID}".sort.cram"
  """
}


// 1.2: MarkDuplicates
process Mrkdup {

  label 'high_mem'
  publishDir path: "${params.outdir}/samples/${sampleID}/picard", mode: "copy", pattern: "*.txt"

  input:
  tuple val(sampleID), val(meta), file(bam), file(bai) from dup_marking

  output:
  file('*.txt') into mrkdup_multiqc
  tuple val(sampleID), val(meta), file('*.md.bam'), file('*.md.bam.bai') into gatk4recaling

  script:
  def taskmem = task.memory == null ? "" : "-Xmx" + javaTaskmem("${task.memory}")
  """
  OUTBAM=\$(echo ${bam} | sed 's/bam/md.bam/')
  OUTMET=\$(echo ${bam} | sed 's/bam/md.metrics.txt/')
  {
  picard ${taskmem} \
    MarkDuplicates \
    TMP_DIR=./ \
    INPUT=${bam} \
    OUTPUT=/dev/stdout \
    COMPRESSION_LEVEL=0 \
    QUIET=TRUE \
    METRICS_FILE=\$OUTMET \
    REMOVE_DUPLICATES=FALSE \
    ASSUME_SORTED=TRUE \
    VALIDATION_STRINGENCY=LENIENT \
    VERBOSITY=ERROR | samtools view -Shb - > \$OUTBAM

  samtools index \$OUTBAM
  } 2>&1 | tee > ${sampleID}.picard_markDuplicates.log.txt
  """
}

// 1.3: GATK4 BestPractices
process Gtkrcl {

  label 'high_mem'
  publishDir path: "${params.outdir}/samples/${sampleID}/gatk4/bestpractice", mode: "copy", pattern: "*.GATK4_BQSR.log.txt "

  input:
  tuple val(sampleID), val(meta), file(bam), file(bai) from gatk4recaling
  file(fasta) from reference.fa
  file(fai) from reference.fai
  file(dict) from reference.dict
  file(dbsnp_files) from reference.dbsnp
  file(intlist) from reference.intlist

  output:
  file('*.table') into gtkrcl_multiqc
  tuple val(sampleID), file('*.bqsr.bam'), file('*.bqsr.bam.bai') into ( gmultimetricing, mosdepthing, mutect2somaticing )
  tuple val(sampleID), val(meta), file('*.bqsr.bam'), file('*.bqsr.bam.bai') into hc_germ
  tuple val(sampleID), val(meta) into metas_pcgr
  file("${sampleID}.GATK4_BQSR.log.txt") into bqsr_log

  script:
  def dbsnp = "${dbsnp_files}/*gz"
  """
  {
  gatk BaseRecalibrator \
    -R ${fasta} \
    -I ${bam} \
    --known-sites \$(echo ${dbsnp}) \
    --use-original-qualities \
    -O ${sampleID}.recal_data.table \
    --disable-sequence-dictionary-validation true \
    -L ${intlist}

  #ApplyBQSR
  OUTBAM=\$(echo ${bam} | sed 's/bam/bqsr.bam/')
  gatk ApplyBQSR \
    -R ${fasta} \
    -I ${bam} \
    --bqsr-recal-file ${sampleID}.recal_data.table \
    --add-output-sam-program-record \
    --use-original-qualities \
    -O \$OUTBAM \
    -L ${intlist}

  samtools index \$OUTBAM \$OUTBAM".bai"
  } 2>&1 | tee >  ${sampleID}.GATK4_BQSR.log.txt
  """
}

// 1.4: scatter-gather implementation for mutect2, octopus
process Scat_gath {

  label 'low_mem'

  input:
  file(intlist) from reference.intlist

  output:
  file('octopus.scatgath.*.bed') into octopus_bedding
  file('mutect2.scatgath.*.bed.interval_list') into mutect2_bedding
  file('hc.scatgath.*.bed.interval_list') into hc_bedding

  script:
  def sgcount = params.scatGath
  if (params.scatGath == null){
    if (params.seqlevel != "wgs"){
      sgcount = 20
    }
    else {
      sgcount = 100
    }
  }
  """
  ##strip out all but chromosomes in the interval_list (no decoys etc)
  CHRS=\$(grep -v "@" ${intlist} | cut -f 1 | uniq)
  for CHR in \$CHRS; do
    grep "SN:\$CHR\\s" ${intlist} >> used.interval_list
  done
  grep -v "@" ${intlist} >> used.interval_list

  ##generate scatters
  picard IntervalListTools \
    I=used.interval_list \
    SCATTER_COUNT=${sgcount} \
    O=./

  ##rename scatters and parse into appropriate format for tools
  ls temp*/* | while read FILE; do
    COUNTN=\$(dirname \$FILE | perl -ane '@s=split(/\\_/); print \$s[1];');
    mv \$FILE mutect2.scatgath.\${COUNTN}.bed.interval_list;
    cp mutect2.scatgath.\${COUNTN}.bed.interval_list hc.scatgath.\${COUNTN}.bed.interval_list
    grep -v @ mutect2.scatgath.\${COUNTN}.bed.interval_list | \
      cut -f 1,2,3,5 > octopus.scatgath.\${COUNTN}.bed
  done
  """
}

process Mosdepth {

  input:
  tuple val(sampleID), file(bam), file(bai) from mosdepthing
  file(bed) from reference.bed

  output:
  file('*') into mosdepth_multiqc

  script:
  """
  mosdepth \
    --no-per-base \
    --by ${bed} \
    ${sampleID} \
    ${bam}
  """
}

/*
================================================================================
                            2.  MUTATION CALLING
================================================================================
*/
// 2.0: GATK4 Germline Haplotypecaller
// Groovy to combine scatter-gather BEDs with bam file for germline
hcbedding = hc_bedding.flatten()
hc_germ
  .map { it -> [it[0],it[1],it[2],it[3]] }
  .combine(hcbedding)
  .set { hcgermbedding }

process Hcaller {

  label 'med_mem'
  errorStrategy 'retry'
  maxRetries 3

  input:
  tuple val(sampleID), val(meta), file(bam), file(bai), file(intlist) from hcgermbedding
  file(fasta) from reference.fa
  file(fai) from reference.fai
  file(dict) from reference.dict
  file(dbsnp_files) from reference.dbsnp
  file(hc_dbs_files) from reference.hc_dbs

  output:
  tuple val(sampleID), val(meta), file('*sort.hc.vcf') into hc_gt

  script:
  def taskmem = task.memory == null ? "" : "--java-options \"-Xmx" + javaTaskmem("${task.memory}") + "\""
  def dbsnp = "${dbsnp_files}/*gz"
  def omni = "${hc_dbs_files}/KG_omni*.gz"
  def kgp1 = "${hc_dbs_files}/KG_phase1*.gz"
  def hpmp = "${hc_dbs_files}/hapmap*.gz"
  """
  SCATGATHN=\$(echo ${intlist} | perl -ane '@s=split(/\\./);print \$s[2];')
  gatk ${taskmem} HaplotypeCaller \
    -R ${fasta} \
    -I ${bam} \
    --dont-use-soft-clipped-bases \
    --standard-min-confidence-threshold-for-calling 20 \
    --dbsnp \$(echo ${dbsnp}) \
    --native-pair-hmm-threads ${task.cpus} \
    -O ${sampleID}".\${SCATGATHN}.hc.vcf" \
    --disable-sequence-dictionary-validation true \
    -L ${intlist}

  picard SortVcf \
    I=${sampleID}".\${SCATGATHN}.hc.vcf" \
    O=${sampleID}".\${SCATGATHN}.sort.hc.vcf" \
    SD=${dict}
  """
}

//group those outputs
hc_gt
  .groupTuple()
  .map { it -> tuple(it[0], it[1].unique(), it[2..-1].flatten()) }
  .set { hc_fm }

// 2.2: HaplotypeCaller merge
process Hc_merge {

  label 'high_mem'
  publishDir path: "${params.outdir}/samples/${sampleID}/haplotypecaller", mode: "copy", pattern: '*.vcf.*'

  input:
  tuple val(sampleID), val(meta), file(rawvcfs) from hc_fm

  output:
  tuple val(sampleID), val(meta), file("${sampleID}.hc.merge.vcf.gz"), file("${sampleID}.hc.merge.vcf.gz.tbi") into ( cpsr_vcf, vep_hc_vcf )

  script:
  """
  ls *.sort.hc.vcf > vcf.list
  picard MergeVcfs I=vcf.list O=${sampleID}".hc.merge.vcf"
  bgzip ${sampleID}".hc.merge.vcf"
  tabix ${sampleID}".hc.merge.vcf.gz"
  """
}

// 2.3: CPSR annotation of GATK4 Germline
process cpsrreport {

  label 'med_mem'

  publishDir "${params.outdir}/reports/cpsr", mode: "copy", pattern: "${metaid}.cpsr.${grchv}.{html,json.gz}"
  publishDir "${params.outdir}/samples/${sampleID}/cpsr", mode: "copy", pattern: "*[!.html]"

  input:
  tuple val(sampleID), val(meta), file(vcf), file(tbi) from cpsr_vcf
  file(grchver) from reference.grchvers
  file(pcgrbase) from reference.pcgrbase

  output:
  file('*') into cpsr_vcfs
  file("${metaid}.cpsr.${grchv}.html") into sendmail_cpsr

  script:
  grchv = "${grchver}".split("\\/")[-1]
  metaid = "${meta}".replaceAll("\\s *", "_").replaceAll("[\\[\\(\\)\\]]","").replaceAll("\"","")
  """
  {
  ##count sample_id as cannot be less than 3
  WCC=\$(echo ${metaid} | tr -d '\\n' | wc -c)
  if (( \$WCC < 3 )); then METAID="_${metaid}";
  else METAID=${metaid}; fi

  ##CPSR v0.6.1
  cpsr.py \
    --no-docker \
    --no_vcf_validate \
    --panel_id 0 \
    --query_vcf ${vcf} \
    --pcgr_dir ${pcgrbase} \
    --output_dir ./ \
    --genome_assembly ${grchv} \
    --conf ${pcgrbase}/data/${grchv}/cpsr_configuration_default.toml \
    --sample_id \$METAID
  } 2>&1 | tee > ${sampleID}.cpsr.log.txt

  for x in \$(ls \$METAID*); do
    nm=\$(echo \$x | sed "s/\$METAID/${metaid}/")
    mv \$x \$nm
  done
  """
}

// 2.4: PicardTools metrics suite for MultiQC HTML report
process Multi_met {

  label 'low_mem'

  publishDir "${params.outdir}/samples/${sampleID}/metrics", mode: "copy"

  input:
  tuple val(sampleID), file(bam), file(bai) from gmultimetricing
  file(fasta) from reference.fa
  file(fai) from reference.fai
  file(dict) from reference.dict
  file(intlist) from reference.intlist

  output:
  file('*.txt') into multimetrics_multiqc

  script:
  def taskmem = task.memory == null ? "" : "-Xmx" + javaTaskmem("${task.memory}")
  """
  {
  if [[ ${params.seqlevel} != "wgs" ]]; then
  picard ${taskmem} CollectHsMetrics \
    I=${bam} \
    O=${sampleID}".hs_metrics.txt" \
    TMP_DIR=./ \
    R=${fasta} \
    BAIT_INTERVALS=${intlist}  \
    TARGET_INTERVALS=${intlist}
  fi
  picard ${taskmem} CollectAlignmentSummaryMetrics \
    I=${bam} \
    O=${sampleID}".AlignmentSummaryMetrics.txt" \
    TMP_DIR=./ \
    R=${fasta}

  picard ${taskmem} CollectMultipleMetrics \
    I=${bam} \
    O=${sampleID}".CollectMultipleMetrics.txt" \
    TMP_DIR=./ \
    R=${fasta}

  picard ${taskmem} CollectSequencingArtifactMetrics \
    I=${bam} \
    O=${sampleID}".artifact_metrics.txt" \
    TMP_DIR=./ \
    R=${fasta}

  picard ${taskmem} CollectInsertSizeMetrics \
    I=${bam} \
    O=${sampleID}".insert_size_metrics.txt" \
    H=${bam}".histogram.pdf" \
    TMP_DIR=./

  } 2>&1 | tee > ${sampleID}.picard.metrics.log
  """
}


mutect2bedding = mutect2_bedding.flatten()
mutect2somaticing
  .map { it -> [it[0],it[1],it[2]]}
  .combine(mutect2bedding)
  .set { mutect2somaticbedding }

// 2.7.1: MuTect2
// NB --germline-resource dollar-sign{dbsnp} removed as no AF causing error

process M2_scatgath {

  label 'med_mem'

  input:
  tuple val(sampleID), file(tumourbam), file(tumourbai), file(intlist) from mutect2somaticbedding
  file(pondir) from reference.hc_dbs
  file(fasta) from reference.fa
  file(fai) from reference.fai
  file(dict) from reference.dict
  file(dbsnp_files) from reference.dbsnp

  output:
  tuple val(sampleID), file('*sort.mutect2.vcf') into mutect2_gt
  tuple val(sampleID), file('*.vcf.stats') into mutect2_st
  tuple val(sampleID), file('*mutect2.f1r2.tar.gz') into mutect2_f1r2
  tuple val(sampleID), file(tumourbam), file(tumourbai) into mutect2_comb

  script:
  def taskmem = task.memory == null ? "" : "--java-options \"-Xmx" + javaTaskmem("${task.memory}") + "\""
  pon = params.assembly == "GRCh37" ? "${pondir}/Mutect2-WGS-panel-b37.vcf.gz" : "${pondir}/KG_pon.hg38.vcf.gz"
  def dbsnp = "${dbsnp_files}/*gz"
  """
  SCATGATHN=\$(echo ${intlist} | perl -ane '@s=split(/\\./);print\$s[2];')
  gatk ${taskmem} \
    Mutect2 \
    --native-pair-hmm-threads ${task.cpus} \
    --reference ${fasta} \
    --input ${tumourbam} \
    --panel-of-normals ${pon} \
    --af-of-alleles-not-in-resource 0.0000025 \
    --output ${sampleID}"."\${SCATGATHN}".mutect2.vcf" \
    --disable-sequence-dictionary-validation true \
    --f1r2-tar-gz \${SCATGATHN}".mutect2.f1r2.tar.gz" \
    -L ${intlist}

  picard SortVcf \
    I=${sampleID}"."\${SCATGATHN}".mutect2.vcf" \
    O=${sampleID}"."\${SCATGATHN}".sort.mutect2.vcf" \
    SD=${dict}
  """
}

// 2.7.2: MuTect2_merge
mutect2_gt
  .groupTuple()
  .map { it -> tuple(it[0], it[1][0..-1].flatten()) }
  .set { mutect2_fm }

process M2_concat {

  label 'med_mem'

  input:
  tuple val(sampleID), file(rawvcfs) from mutect2_fm

  output:
  tuple val(sampleID), file('*mutect2.merge.vcf') into mutect2_merge

  script:
  """
  ls *.sort.mutect2.vcf > vcf.list
  picard MergeVcfs I=vcf.list O=${sampleID}".mutect2.merge.vcf"
  """
}

mutect2_st
  .groupTuple()
  .map { it -> tuple(it[0], it[1][0..-1].flatten()) }
  .set { mutect2_sm }

// 2.7.3: MuTect2 Concatenate VCFs
process M2_conc_stat {

  label 'med_mem'

  input:
  tuple val(sampleID), file(stats) from mutect2_sm

  output:
  tuple val(sampleID), file('*mutect2.merge.vcf.stats') into mutect2_stats

  script:
  """
  STATS=\$(ls *stats | perl -ane 'foreach \$k (@F){print "--stats \$k ";}')
  gatk MergeMutectStats --output ${sampleID}".mutect2.merge.vcf.stats" \$STATS
  """
}

// 2.7.4: MuTect2 Concatenate VCFs
mutect2_f1r2.groupTuple()
            .map { it -> [it[0], it[1..-1].flatten()] }
            .set { mutect2_f1r2_set }

process M2_f1r2_comb {

  label 'med_mem'

  input:
  tuple val(sampleID), file(mutect2_ro) from mutect2_f1r2_set

  output:
  tuple val(sampleID), file("${sampleID}.mutect2.f1r2.tar.gz") into mutect2_f1r2_comb

  script:
  """
  ALL_F1R2_INPUT=\$(for x in *.mutect2.f1r2.tar.gz; do echo -n "-I \$x "; done)
  gatk LearnReadOrientationModel \$ALL_F1R2_INPUT -O ${sampleID}.mutect2.f1r2.tar.gz
  """
}

// 2.7.5: MuTect2 Contamination
mutect2_comb
  .join(mutect2_merge)
  .join(mutect2_stats)
  .join(mutect2_f1r2_comb)
  .groupTuple()
  .map { it -> [it[0], it[1..5].flatten()].flatten() }
  .set { mutect2_stats_merge }

process M2_con_filt {

  label 'med_mem'

  publishDir path: "${params.outdir}/samples/${sampleID}/mutect2", mode: "copy", overwrite: true

  input:
  tuple val(sampleID), file(tumourbam), file(tumourbai), file(mergevcf), file(statsvcf), file(readorient) from mutect2_stats_merge
  file(fasta) from reference.fa
  file(fai) from reference.fai
  file(dict) from reference.dict
  file(levelbase) from reference.seqlevel
  file(intlist) from reference.intlist

  output:
  tuple val(sampleID), file('*snv_indel.pass.vcf') into mutect2_veping
  file('*') into completedmutect2call

  script:
  def taskmem = task.memory == null ? "" : "--java-options \"-Xmx" + javaTaskmem("${task.memory}") + "\""
  hg = params.assembly == "GRCh37" ? "hg19" : "hg38"
  gpsgz = params.seqlevel != "wgs" ? "${levelbase}/${params.levelTag}/af-only-gnomad.${params.levelTag}.${hg}.noChr.vcf.gz" : "${levelbase}/af-only-gnomad.wgs.${hg}.noChr.vcf.gz"
  """
  gatk ${taskmem} \
    GetPileupSummaries \
    -I ${tumourbam} \
    -V ${gpsgz} \
    -O ${sampleID}".getpileupsummaries.table" \
    -L ${intlist}

  gatk CalculateContamination \
    -I ${sampleID}".getpileupsummaries.table" \
    -O ${sampleID}".calculatecontamination.table"

  CONTAM=\$(tail -n+2 ${sampleID}.calculatecontamination.table | cut -f 2 | cut -d "." -f 1)
  if [[ \$CONTAM != 0 ]]; then
    touch ${sampleID}".CONTAMINATION.WARNING.txt"
  fi

  gatk IndexFeatureFile \
    --input ${mergevcf}

  gatk ${taskmem} \
    FilterMutectCalls \
    --reference ${fasta} \
    --contamination-table ${sampleID}".calculatecontamination.table" \
    --interval-padding 5 \
    --output ${sampleID}".mutect2.FilterMutectCalls.vcf" \
    --unique-alt-read-count 3 \
    --variant ${mergevcf} \
    --stats ${statsvcf} \
    --disable-sequence-dictionary-validation true \
    --ob-priors ${readorient} \
    -L ${intlist}

  perl ${workflow.projectDir}/bin/filter_Lancet_Mutect2_Manta-Strelka2_Format.pl \
    ID=${sampleID} \
    DP=14 \
    MD=2 \
    VCF=${sampleID}".mutect2.FilterMutectCalls.vcf"
  """
}

/*
================================================================================
                          3.  ANNOTATION AND REPORTING
================================================================================
*/

// 3.01: HC_merge VEP
process vepHC {

  label 'low_mem'

  publishDir path: "${params.outdir}/samples/${sampleID}/haplotypecaller", mode: "copy"

  input:
  tuple val(sampleID), val(meta), file(vcf), file(tbi) from vep_hc_vcf
  file(fasta) from reference.fa
  file(fai) from reference.fai
  file(dict) from reference.dict
  file(grchver) from reference.grchvers
  file(pcgrbase) from reference.pcgrbase

  output:
  tuple val(sampleID), val(meta), file("${sampleID}.hc.merge.vep.vcf") into hc_vepd

  script:
  grch_vers = "${grchver}".split("\\/")[-1]
  """
  vep --dir_cache ${pcgrbase}/data/${grch_vers}/.vep \
    --offline \
    --assembly ${params.assembly} \
    --vcf_info_field ANN \
    --symbol \
    --species homo_sapiens \
    --check_existing \
    --cache \
    --fork ${task.cpus} \
    --af_1kg \
    --af_gnomad \
    --vcf \
    --input_file ${vcf} \
    --output_file ${sampleID}.hc.merge.vep.vcf \
    --format "vcf" \
    --fasta ${fasta} \
    --hgvs \
    --canonical \
    --ccds \
    --force_overwrite \
    --verbose
  """
}

process vepann {

  label 'med_mem'

  publishDir path: "${params.outdir}/samples/${sampleID}/mutect2", mode: "copy", pattern: "${sampleID}.mutect2.snv_indel.pass.vep.vcf"

  input:
  tuple val(sampleID), file(vcf) from mutect2_veping
  file(fasta) from reference.fa
  file(fai) from reference.fai
  file(dict) from reference.dict
  file(grchver) from reference.grchvers
  file(pcgrbase) from reference.pcgrbase

  output:
  tuple val(sampleID), file("*.vep.vcf") into runPCGR

  script:
  def grch_vers = "${grchver}".split("\\/")[-1]
  def vcf_anno = "${vcf}".replaceAll(".vcf", ".vep.vcf")
  """
  vep --dir_cache ${pcgrbase}/data/${grch_vers}/.vep \
    --offline \
    --assembly ${params.assembly} \
    --vcf_info_field ANN \
    --symbol \
    --species homo_sapiens \
    --check_existing \
    --cache \
    --fork ${task.cpus} \
    --af_1kg \
    --af_gnomad \
    --vcf \
    --input_file ${vcf} \
    --output_file ${vcf_anno} \
    --format "vcf" \
    --fasta ${fasta} \
    --hgvs \
    --canonical \
    --ccds \
    --force_overwrite \
    --verbose

  ##remove calls on decoys
  perl -ane 'if(\$F[0]=~m/^#/){print \$_;}else{if(\$F[0]!~m/_/){print \$_;}}' ${vcf_anno} >1
  mv 1 ${vcf_anno}
  """
}

// 3.3 PCGR report
// take all mutations in consensus.tab from pass.vcfs into single VCF for PCGR

runPCGR
  .join(pcgr_meta)
  .groupTuple()
  .map { it -> [it[0], it[1], it[2]].flatten() }
  .set { run_pcgr_meta }

process pcgrreport {

  label 'low_mem'
  errorStrategy 'retry'
  maxRetries 3

  publishDir "${params.outdir}/reports/pcgr", mode: "copy", pattern: "${sampleID}.pcgr_acmg.${grch_vers}.{html,json.gz}"
  publishDir "${params.outdir}/samples/${sampleID}/pcgr", mode: "copy"

  input:
  tuple val(sampleID), file(vcf), val(meta) from run_pcgr_meta
  file(grchver) from reference.grchvers
  file(pcgrbase) from reference.pcgrbase
  file(levelbase) from reference.seqlevel

  output:
  file('*') into completed_pcgr
  tuple file("${metaid}.pcgr_acmg.${grch_vers}.html"), file("${metaid}.pcgr_acmg.${grch_vers}.json.gz") into sendmail_pcgr

  script:
  grch_vers = "${grchver}".split("\\/")[-1]
  config = params.seqlevel != "wgs" ? "${levelbase}/${params.levelTag}/pcgr_configuration_${params.levelTag}.toml" : "${pcgrbase}/data/${grch_vers}/pcgr_configuration_default.toml"
  metaid = "${meta}".replaceAll("\\s *", "_").replaceAll("[\\[\\(\\)\\]]","").replaceAll("\"","")
  assay = params.seqlevel == "wgs" ? "WGS" : params.seqlevel == "exome" ? "WES" : "TARGETED"
  tmb_msi = params.seqlevel == "panel" ? "" : "--estimate_tmb --estimate_msi_status --tmb_algorithm all_coding"
  """
  {
  ##count sample_id as cannot be less than 3
  WCC=\$(echo ${metaid} | tr -d '\\n' | wc -c)
  if (( \$WCC < 3 )); then METAID="_${metaid}";
  else METAID=${metaid}; fi

  ##PCGR 0.9.1
  pcgr.py \
    --pcgr_dir ${pcgrbase} \
    --output_dir ./ \
    --genome_assembly ${grch_vers} \
    --conf ${config} \
    --sample_id \$METAID \
    --input_vcf ${vcf} \
    --no-docker \
    --force_overwrite \
    --no_vcf_validate \
    --tumor_only \
    --include_trials \
    --assay ${assay} ${tmb_msi}

  for x in \$(ls \$METAID*); do
    nm=\$(echo \$x | sed "s/\$METAID/${metaid}/")
    mv \$x \$nm
  done

  } 2>&1 | tee > ${sampleID}.pcgr.log.txt
  """
}

/*
================================================================================
                          4.  MULTIQC AND CLOSEOUT
================================================================================
*/
// 4.0 Run multiQC to finalise report
process MultiQC {

  label 'low_mem'
  publishDir path: "${params.outdir}/reports/multiQC", mode: "copy"

  input:
  file(fastps) from fastp_multiqc.collect()
  file(fastqcs) from fastqc_multiqc.collect()
  file(gtkrcls) from gtkrcl_multiqc.collect()
  file(multimetrics) from multimetrics_multiqc.collect()
  file(mrkdups) from mrkdup_multiqc.collect()
  file(mosdepth) from mosdepth_multiqc.collect()

  output:
  file('*') into completedmultiqc
  file("*.html") into sendmail_multiqc

  script:
  """
  multiqc . -i ${params.runID}".tumour_only" --tag DNA -f -c ${params.multiqcConfig}
  """
}

// 4.1.1: container software versions
process software_vers {

  label 'low_mem'
  publishDir "pipeline_info", mode: 'copy'
  publishDir path: "${params.outdir}/reports/software_vers", mode: "copy"

  output:
  file('*') into ch_software_vers

  script:
  """
  conda env export > somenone_software_versions.yaml
  pcgr.py --version > pcgr_software_versions.txt
  cpsr.py --version >> pcgr_software_versions.txt
  """
}

// 4.2: ZIP for sending on sendmail

sendmail_pcgr
  .mix(sendmail_cpsr)
  .mix(sendmail_multiqc)
  .set { sendmail_all }

process zipup {

  label 'low_mem'

  input:
  file(send_all) from sendmail_all.collect()

  output:1
  file("${params.runID}.tumour_only.zip") into send_zip

  script:
  """
  mkdir -p ${params.runID}/html_reports && mv *html ${params.runID}/html_reports/
  if [[ \$(find ./ -maxdepth 0 | wc -l) > 1 ]]; then
    mkdir other && mv *.* ./other/
  fi
  zip -r ${params.runID}.tumour_only.zip ${params.runID}/html_reports
  """
}

// 4.3: Completion e-mail notification

workflow.onComplete {
  sleep(1000)
  def subject = """\
    [brucemoran/tumour_only] SUCCESS: $params.runID [$workflow.runName]
    """
    .stripIndent()
  if (!workflow.success) {
      subject = """\
        [brucemoran/tumour_only] FAILURE: $params.runID [$workflow.runName]
        """
        .stripIndent()
  }

  def msg = """\
    Pipeline execution summary
    ---------------------------
    RunID       : ${params.runID}
    RunName     : ${workflow.runName}
    Completed at: ${workflow.complete}
    Duration    : ${workflow.duration}
    workDir     : ${workflow.workDir}
    exit status : ${workflow.exitStatus}
    """
    .stripIndent()

  def attachments = send_zip.toList().getVal()

  sendMail(to: "${params.email}",
           subject: subject,
           body: msg,
           attach: attachments)
}
