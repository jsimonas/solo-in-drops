#!/usr/bin/env nextflow
/*
========================================================================================
                                    solo-in-drops
========================================================================================
 jsimonas/solo-in-drops Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/jsimonas/solo-in-drops
----------------------------------------------------------------------------------------
*/

def helpMessage() {
    // TODO nf-core: Add to this help message with new command line parameters
    log.info nfcoreHeader()
    log.info"""

    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run jsimonas/solo-in-drops --run_dir 'path/to/bcl_folder' --sample_sheet 'path/to/extended_sample_sheet.xlsx' -profile singularity

    Mandatory arguments:
      --run_dir [path/to/folder]      Path to input data (must be surrounded with quotes)
      --run_module [str]              Pipeline module to run. Can be set as "complete", "demux" or "fastq". If the latter selected, sample sheet is not required. Default: "complete".
      --scrna_protocol [str]          Protocol used to generate scRNA-seq libraries. Default: "indrops". Currently, customized "indrops" or "splitpool" protocols are supported.
      --sample_sheet [file]           Full path to extended sample sheet file. Example can be found at solo-in-drops/assets/extended_sample_sheet_template.xlsx
      --sequencer [str]               Sequencer used to generate the data. Default: "nextseq". Can be set as "nextseq" or "miseq".
      --align_mode [str]              STAR alignment mode. Default: "cell". Can be set as "bacteria" to switch off splice alignments.
      -profile [str]                  Configuration profile to use. Can use multiple (comma separated).
                                      Available: conda, docker and singularity

    References:                       If not specified in the configuration file or you wish to overwrite any of the references
      --star_index [path/to/folder]   Path to star index directory (same as --genomeDir parameter in STAR)
      --barcode_whitelist [file]      Path to cell barcode list (a text file containing one barcode sequence per line)
    
    STARsolo arguments:               If not specified, the default parameters will be used
      --bc_read_length [int]          Read length of cell barcode read. Default: equal to sum of BC + UMI
      --solo_multi_mappers            Allow multi-gene read quantification. Can be set as "Uniform", "PropUnique", "EM", "Rescue" or any combination of these options. Default: "Uniform"
      --solo_features                 Counting option for transcriptomic features, including "Gene" (default), "GeneFull", and "Velocyto" or any combination of these options. Combinations should be provided in one string, i.e. "Gene Velocyto".           

    bcl2fastq arguments:              If not specified, the default parameters will be used
      --barcode_mismatches            Allowed missmaches in library index for demultiplexing. Default: 1
      --write_fastq                   Output raw FASTQ files for each library. Default: false
    
    Other options:
      --outdir [file]                 The output directory where the results will be saved
      -name [str]                     Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic
    """.stripIndent()
}

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

/*
 * SET UP CONFIGURATION VARIABLES
 */

// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
if (!(workflow.runName ==~ /[a-z]+_[a-z]+/)) {
    custom_runName = workflow.runName
}

// Validate mandatory inputs

if (!(params.run_module.equals('complete') || params.run_module.equals('demux') || params.run_module.equals('fastq'))){
    exit 1, "Uncorrect pipeline run module was provided! Can be set as 'complete', 'demux' or 'fastq' module."
}

if (params.sample_sheet && (params.run_module.equals('complete') || params.run_module.equals('demux'))){
    sheet_file = file(params.sample_sheet, checkIfExists: true)
} else if (params.run_module.equals('fastq')) {
    sheet_file = Channel.empty()
} else {
    exit 1, "The extended sample sheet is not provided! Template of the file can be found at solo-in-drops/assets/extended_sample_sheet_template.xlsx"
}

// Check run directory
if (params.run_dir){
    runDir = file(params.run_dir, checkIfExists: true)
    } else {
    exit 1, "Input directory not found!"
    }
runName = runDir.getName()

if (!(params.sequencer.equals('nextseq') || params.sequencer.equals('miseq'))){
    exit 1, "Unsupported sequencer provided! Currently nextseq or miseq are supported"
}

// Check STAR index
if( params.star_index ){
    star_index = Channel
        .fromPath(params.star_index)
        .ifEmpty { exit 1, "STAR index not found: ${params.star_index}" }
}

//Check barcode whitelist
if( params.barcode_whitelist ){
    barcode_whitelist = Channel
        .fromFilePairs(params.barcode_whitelist)
        .ifEmpty { exit 1, "barcode whitelist not found: ${params.barcode_whitelist}" }
}
if( params.barcode_whitelist ){
    bbarcode_whitelist = Channel
        .fromFilePairs(params.barcode_whitelist)
        .ifEmpty { exit 1, "barcode whitelist not found: ${params.barcode_whitelist}" }
}
bbarcode_whitelist.flatten().view()

// bc_wl = whitelist.join(' ')


// Define scRNA protocol related parameters

// bcl2fastq
if (!params.scrna_protocol.equals("splitpool")){
    mask = 'y*,I*,y*,y*'
} else {
    mask = 'y*,I*,y*'
}
// STARsolo
if (!params.scrna_protocol.equals("splitpool")){
    cb_position = '0_0_0_7 0_8_0_15'
    umi_position = '0_16_0_23'
} else {
    cb_position = '0_0_0_7 0_8_0_17 0_18_0_25'
    cb_position = '0_26_0_33'
}

if (!(params.align_mode.equals('cell') || params.align_mode.equals('bacteria'))){
    exit 1, "Provided alingment mode is not supported! Please use 'cell' or 'bacteria' for --align_mode parameter."
} else if (params.align_mode.equals('cell')){
    twopass = 'Basic'
    alignintronmax = 0
} else if (params.align_mode.equals('bacteria')){
    twopass = 'None'
    alignintronmax = 1
}

// Stage config files
ch_multiqc_config = file("$baseDir/assets/multiqc_config.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
ch_output_docs = file("$baseDir/docs/output.md", checkIfExists: true)

// Header log info
log.info nfcoreHeader()
def summary = [:]
if (workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Run Name']         = custom_runName ?: workflow.runName
// TODO nf-core: Report custom parameters here
summary['Run module']       = params.run_module
summary['Sample sheet']     = params.sample_sheet
summary['Input directory']  = params.run_dir
summary['Sequencer']        = params.sequencer
summary['scRNAseq protocol']     = params.scrna_protocol
summary['Aligment mode']    = params.align_mode
summary['Multi-mapper recovery'] = params.solo_multi_mappers
summary['STARsolo features'] = params.solo_features
summary['STAR index']       = params.star_index
summary['CB whitelist']     = params.barcode_whitelist
summary['Max Resources']    = "$params.max_memory memory, $params.max_cpus cpus, $params.max_time time per job"
if (workflow.containerEngine) summary['Container'] = "$workflow.containerEngine - $workflow.container"
summary['Output dir']       = params.outdir
summary['Launch dir']       = workflow.launchDir
summary['Working dir']      = workflow.workDir
summary['Script dir']       = workflow.projectDir
summary['User']             = workflow.userName
summary['Config Profile'] = workflow.profile
if (params.config_profile_description) summary['Config Description'] = params.config_profile_description
if (params.config_profile_contact)     summary['Config Contact']     = params.config_profile_contact
if (params.config_profile_url)         summary['Config URL']         = params.config_profile_url
if (params.email || params.email_on_fail) {
    summary['E-mail Address']    = params.email
    summary['E-mail on failure'] = params.email_on_fail
    summary['MultiQC maxsize']   = params.max_multiqc_email_size
}
log.info summary.collect { k,v -> "${k.padRight(18)}: $v" }.join("\n")
log.info "-\033[2m--------------------------------------------------\033[0m-"

// Check the hostnames against configured profiles
checkHostname()

Channel.from(summary.collect{ [it.key, it.value] })
    .map { k,v -> "<dt>$k</dt><dd><samp>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }
    .reduce { a, b -> return [a, b].join("\n            ") }
    .map { x -> """
    id: 'solo-in-drops-summary'
    description: " - this information is collected when the pipeline is started."
    section_name: 'jsimonas/solo-in-drops Workflow Summary'
    section_href: 'https://github.com/jsimonas/solo-in-drops'
    plot_type: 'html'
    data: |
        <dl class=\"dl-horizontal\">
            $x
        </dl>
    """.stripIndent() }
    .set { ch_workflow_summary }

/*
 * Parse software version numbers
 */
process get_software_versions {
    publishDir "${params.outdir}/pipeline_info", mode: 'copy',
//        saveAs: { filename -> "$runName"+"_"+"$filename" }      
        saveAs: { filename ->
                      if (filename.indexOf(".csv") > 0) "$runName"+"_"+"$filename" 
                      else null
                }

    output:
    file 'software_versions_mqc.yaml' into ch_software_versions_yaml
    file "software_versions.csv"

    script:
    // versions of tools
    """
    echo $workflow.manifest.version > v_pipeline.txt
    echo $workflow.nextflow.version > v_nextflow.txt
    STAR --version &> v_star.txt 2>&1
    samtools --version |& grep "sam" &> v_samtools.txt
    bcl2fastq --version |& grep "bcl" &> v_bcl2fastq.txt
    seqkit version &> v_seqkit.txt
    fastqc -version |& grep "v" &> v_fastqc.txt
    multiqc --version &> v_multiqc.txt || true
    scrape_software_versions.py &> software_versions_mqc.yaml
    """
}

/*
 * STEP 1 - convert extended to standard sample sheet
 */
process convert_sample_sheet {
    tag "$sheet"
    label 'process_low'
    publishDir path: "${params.outdir}/", mode: 'copy'
 
    input:
    file sheet from sheet_file

    when:
    params.run_module.equals('complete') || params.run_module.equals('demux') 

    output:
    file "*.csv" into standard_samplesheet
    
    script:
    """
    convert_to_samplesheet.py --file "${sheet}" --out "standard_samplesheet.csv"
    """
}

/*
 * STEP 2 - convert bcl to fastq files
 */
process bcl_to_fastq {
    tag "$runName"
    label 'process_high'
    publishDir path: "${params.outdir}/", pattern: "*/**.fastq.gz", mode: 'copy',
        saveAs: { filename -> 
            if (params.write_fastq) filename
            else null
        }
    publishDir path: "${params.outdir}/", pattern: "*.fastq.gz", mode: 'copy',
        saveAs: { filename -> 
            if (params.write_fastq) "Undetermined/fastqs/$filename"
            else null
        }

    input:
    file sheet from standard_samplesheet

    when:
    params.run_module.equals('complete') || params.run_module.equals('demux') 

    output:
    file "*/**{R1,R2,R3}_001.fastq.gz" into fastqs_fqc_ch, fastqs_output_ch mode flatten
    file "*{R1,R2,R3}_001.fastq.gz" into und_fastqs_fqc_ch mode flatten
    file "Stats" into bcl2fq_stats_ch mode flatten

    script:
    
    """
    bcl2fastq \\
    --runfolder-dir ${runDir} \\
    --output-dir . \\
    --sample-sheet ${sheet} \\
    --mask-short-adapter-reads 0 \\
    --minimum-trimmed-read-length 0 \\
    --use-bases-mask ${mask} \\
    --no-lane-splitting \\
    --create-fastq-for-index-reads \\
    --barcode-mismatches $params.barcode_mismatches \\
    --processing-threads $task.cpus
    """
}

// add project name
fqname_fqfile_ch = fastqs_fqc_ch.map{
    file -> [file.getParent().getName(), file]
}
undetermined_fqfile_ch = und_fastqs_fqc_ch.map{
    file -> ["Undetermined", file]
}
fastqs = Channel.empty()
fastqs_ch = fastqs.mix(fqname_fqfile_ch, undetermined_fqfile_ch)

/*
 * STEP 3 - FastQC
 */
process fastqc {
    tag "$projectName"
    label 'process_medium'
    publishDir "${params.outdir}/${projectName}/fastqc", mode: 'copy'
        
    input:
    set val(projectName), file(fastq) from fastqs_ch

    when:
    params.run_module.equals('complete') || params.run_module.equals('demux')

    output:
    set val(projectName), file("*_fastqc.{zip,html}") into fastqc_results_ch

    script:
    """
    fastqc --quiet --threads $task.cpus ${fastq}
    """
}

// make paired channel for fastqs
fastqs_output_ch.flatMap()
            .map{ file ->
                if ( "${file}".contains("_R1_") || "${file}".contains("_R2_") || "${file}".contains("_R3_")){
                    def key_match = file.name.toString() =~ /(.+)_R\d+_001\.fastq\.gz/
                    def key = key_match[0][1]
                    def proj = file.getParent().getName()
                    return tuple(key, proj, file)
                }
            }
            .groupTuple(by: [0,1])
            .set{ fastq_pairs_ch }

/*
 * STEP 4 - Merge FASTQ
 */
process mergefastq {
    tag "$prefix"
    label 'process_high'
    publishDir "${params.outdir}/${projectName}/merged_fastq", mode: 'copy'
    echo true
    
    input:
    set val(prefix), val(projectName), file(reads) from fastq_pairs_ch
    
    when:
    params.run_module.equals('complete') || params.run_module.equals('demux') 

    output:
    set val(prefix), val(projectName), file('*_{bc,cdna}_001.fastq.gz') into merged_fastq_ch
    
    script:
    R1 = reads[0]
    R2 = reads[1]
    R3 = reads[2]
    
    if (params.scrna_protocol.equals('indrops') && params.sequencer.equals('miseq')){
    """
    seqkit concat ${R2} ${R1} \\
    --out-file ${prefix}_bc_001.fastq.gz \\
    --line-width 0 \\
    --threads $task.cpus
    cp ${R3} ${prefix}_cdna_001.fastq.gz
    """
    } else if (params.scrna_protocol.equals('indrops') && params.sequencer.equals('nextseq')){
    """
    seqkit concat <(seqkit seq --reverse --complement --seq-type 'dna' ${R2}) ${R1} \\
    --out-file ${prefix}_bc_001.fastq.gz \\
    --line-width 0 \\
    --threads $task.cpus
    cp ${R3} ${prefix}_cdna_001.fastq.gz
    """
    } else if (params.scrna_protocol.equals('splitpool')){
    """
    zcat ${R1} \\
    | awk 'NR%4==2 || NR%4==0{\$0=substr(\$0,5,8)substr(\$0,18,10)substr(\$0,32,8)substr(\$0,1,4)substr(\$0,40,4)}1 ' \\
    | gzip > ${prefix}_bc_001.fastq.gz
    cp ${R2} ${prefix}_cdna_001.fastq.gz
    """
    }
}

// assign merged fastq channel
if(params.run_module.equals('fastq')){
    merged_fastq_paired_ch = Channel
        .fromFilePairs("$runDir/*_{bc,cdna}_001.fastq.gz", size: -1)
        { file -> tags = (file.name =~ /(.+)(_S\d+)_\S+_001/)[0]; tags[1]+tags[2]+","+tags[1] }
        .ifEmpty {
            error "Cannot find any reads matching bc_001.fastq.gz and cdna_001.fastq.gz in the: ${params.run_dir}"
        }
        .map {
            tag, pair -> tags = tag.split(/,/) ; [tags[0], tags[1], pair] 
        }
} else {
    merged_fastq_paired_ch = merged_fastq_ch
}

/*
 * STEP 5 - STARsolo
 */
process starsolo {
    tag "$prefix"
    label 'process_high'
    publishDir "${params.outdir}/", mode: 'copy',
        saveAs: {
            filename -> 
            if(params.run_module.equals('fastq')){
                "starsolo/$prefix/$filename"
            }
            else {
                "${projectName}/starsolo/$prefix/$filename"
            }
        }
    echo true

    input:
    set val(prefix), val(projectName), file(reads) from merged_fastq_paired_ch
    file(whitelist) from barcode_whitelist
    file index from star_index.collect()

    when:
    !(params.run_module.equals('demux')) 

    output:
    file "*.bam"
    file "*.out" 
    set val(projectName), file("*.final.out") into alignment_logs
    set val(projectName), file("*_Solo.out/*/${prefix}_*.{stats,txt,csv}") into features_stats_ch
    set val(projectName), file("*_Solo.out/${prefix}_Barcodes.stats") into barcodes_stats_ch

    script:
    prefix = reads[0].toString() - ~/(_bc_001)?(\.fastq)?(\.gz)?$/
    bc_read = reads[0]
    cdna_read = reads[1]
    bc_wl = whitelist.join(' ')
    
    """
    STAR \\
    --genomeDir ${index} \\
    --readFilesIn ${cdna_read} ${bc_read} \\
    --soloCBwhitelist ${bc_wl} \\
    --runThreadN ${task.cpus} \\
    --outFileNamePrefix ${prefix}_ \\
    --alignIntronMax ${alignintronmax} \\
    --outSAMunmapped Within \\
    --outSAMtype BAM SortedByCoordinate \\
    --outBAMsortingBinsN 20 \\
    --outSAMattributes NH HI nM AS CR UR CB UB sS sQ sM GX GN \\
    --twopassMode ${twopass} \\
    --runDirPerm All_RWX \\
    --readFilesCommand zcat \\
    --soloMultiMappers ${params.solo_multi_mappers} \\
    --soloFeatures ${params.solo_features} \\
    --soloType CB_UMI_Complex \\
    --soloCBposition ${cb_position} \\
    --soloUMIposition ${umi_position} \\
    --soloBarcodeReadLength ${params.bc_read_length} \\
    --soloCBmatchWLtype EditDist_2 
    
    feature="\$(echo "$params.solo_features" | sed 's/\s.*\$//')"
    
    awk '{print NR "\t" \$0}' ${prefix}_Solo.out/\${feature}/UMIperCellSorted.txt \\
    > ${prefix}_Solo.out/\${feature}/${prefix}_UMIperCellSorted.txt
    
    awk 'gsub(/^\s+/,"", \$0)gsub(/\s+/,"\t")' ${prefix}_Solo.out/\${feature}/Features.stats \\
    > ${prefix}_Solo.out/\${feature}/${prefix}_Features.stats
    
    awk 'gsub(/^\s+/,"", \$0)gsub(/\s+/,"\t")' ${prefix}_Solo.out/Barcodes.stats \\
    > ${prefix}_Solo.out/${prefix}_Barcodes.stats
    
    awk 'BEGIN{FS=","}{for (i=1;i<=NF;i++) col[i] = col[i]","\$i} END{for (i=1;i<=NF;i++) print col[i]}' ${prefix}_Solo.out/\${feature}/Summary.csv \\
    | awk 'BEGIN{FS=OFS=""}NR==1{print "Sample Name" OFS \$0}NR==2{print "${prefix}" OFS \$0}' \\
    >  ${prefix}_Solo.out/"\${feature}"/${prefix}_Summary.csv
    
    """
}


/*
 * STEP 6 - MultiQC 
 */
process multiqc {
    tag "$runName"
    label 'process_low'
    publishDir "${params.outdir}/multiqc", mode: 'copy',
      saveAs: { filename -> "${runName}_${filename}" }

    input:
    file (multiqc_config) from ch_multiqc_config
    file (mqc_custom_config) from ch_multiqc_custom_config.collect().ifEmpty([])
    file bcl2fq_stats from bcl2fq_stats_ch.collect().ifEmpty([])
    file (fastqc:"fastqc/*") from fastqc_results_ch.flatten().collect().ifEmpty([])
    file (starsolo:"starsolo/*") from alignment_logs.collect().ifEmpty([])
    file (starsolo:"starsolo/*") from features_stats_ch.flatten().collect().ifEmpty([])
    file (starsolo:"starsolo/*") from barcodes_stats_ch.collect().ifEmpty([])
    file workflow_summary from ch_workflow_summary.collectFile(name: "workflow_summary_mqc.yaml")
    file ("software_versions/*") from ch_software_versions_yaml.collect()
    
    output:
    file "*multiqc_report.html" into ch_multiqc_report
    file "*_data"
    file "multiqc_plots"

    script:
    rtitle = custom_runName ? "--title \"$custom_runName\"" : ''
    rfilename = custom_runName ? "--filename " + custom_runName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report" : ''
    custom_config_file = params.multiqc_config ? "--config $mqc_custom_config" : ''
    """
    multiqc -f $rtitle $rfilename $custom_config_file .
    """
}

/*
 * STEP 7 - Output Description HTML
 */
process output_documentation {
    publishDir "${params.outdir}/pipeline_info", mode: 'copy',
        saveAs: { filename -> "${runName}_${filename}" }

    input:
    file output_docs from ch_output_docs

    output:
    file "results_description.html"

    script:
    """
    markdown_to_html.py $output_docs -o results_description.html
    """
}

/*
 * Completion e-mail notification
 */
workflow.onComplete {

    // Set up the e-mail variables
    def subject = "[jsimonas/solo-in-drops] Successful: $workflow.runName"
    if (!workflow.success) {
        subject = "[jsimonas/solo-in-drops] FAILED: $workflow.runName"
    }
    def email_fields = [:]
    email_fields['version'] = workflow.manifest.version
    email_fields['runName'] = custom_runName ?: workflow.runName
    email_fields['success'] = workflow.success
    email_fields['dateComplete'] = workflow.complete
    email_fields['duration'] = workflow.duration
    email_fields['exitStatus'] = workflow.exitStatus
    email_fields['errorMessage'] = (workflow.errorMessage ?: 'None')
    email_fields['errorReport'] = (workflow.errorReport ?: 'None')
    email_fields['commandLine'] = workflow.commandLine
    email_fields['projectDir'] = workflow.projectDir
    email_fields['summary'] = summary
    email_fields['summary']['Date Started'] = workflow.start
    email_fields['summary']['Date Completed'] = workflow.complete
    email_fields['summary']['Pipeline script file path'] = workflow.scriptFile
    email_fields['summary']['Pipeline script hash ID'] = workflow.scriptId
    if (workflow.repository) email_fields['summary']['Pipeline repository Git URL'] = workflow.repository
    if (workflow.commitId) email_fields['summary']['Pipeline repository Git Commit'] = workflow.commitId
    if (workflow.revision) email_fields['summary']['Pipeline Git branch/tag'] = workflow.revision
    email_fields['summary']['Nextflow Version'] = workflow.nextflow.version
    email_fields['summary']['Nextflow Build'] = workflow.nextflow.build
    email_fields['summary']['Nextflow Compile Timestamp'] = workflow.nextflow.timestamp

    // TODO nf-core: If not using MultiQC, strip out this code (including params.max_multiqc_email_size)
    // On success try attach the multiqc report
    def mqc_report = null
    try {
        if (workflow.success) {
            mqc_report = ch_multiqc_report.getVal()
            if (mqc_report.getClass() == ArrayList) {
                log.warn "[jsimonas/solo-in-drops] Found multiple reports from process 'multiqc', will use only one"
                mqc_report = mqc_report[0]
            }
        }
    } catch (all) {
        log.warn "[jsimonas/solo-in-drops] Could not attach MultiQC report to summary email"
    }

    // Check if we are only sending emails on failure
    email_address = params.email
    if (!params.email && params.email_on_fail && !workflow.success) {
        email_address = params.email_on_fail
    }

    // Render the TXT template
    def engine = new groovy.text.GStringTemplateEngine()
    def tf = new File("$baseDir/assets/email_template.txt")
    def txt_template = engine.createTemplate(tf).make(email_fields)
    def email_txt = txt_template.toString()

    // Render the HTML template
    def hf = new File("$baseDir/assets/email_template.html")
    def html_template = engine.createTemplate(hf).make(email_fields)
    def email_html = html_template.toString()

    // Render the sendmail template
    def smail_fields = [ email: email_address, subject: subject, email_txt: email_txt, email_html: email_html, baseDir: "$baseDir", mqcFile: mqc_report, mqcMaxSize: params.max_multiqc_email_size.toBytes() ]
    def sf = new File("$baseDir/assets/sendmail_template.txt")
    def sendmail_template = engine.createTemplate(sf).make(smail_fields)
    def sendmail_html = sendmail_template.toString()

    // Send the HTML e-mail
    if (email_address) {
        try {
            if (params.plaintext_email) { throw GroovyException('Send plaintext e-mail, not HTML') }
            // Try to send HTML e-mail using sendmail
            [ 'sendmail', '-t' ].execute() << sendmail_html
            log.info "[solo-in-drops] Sent summary e-mail to $email_address (sendmail)"
        } catch (all) {
            // Catch failures and try with plaintext
            [ 'mail', '-s', subject, email_address ].execute() << email_txt
            log.info "[solo-in-drops] Sent summary e-mail to $email_address (mail)"
        }
    }

    // Write summary e-mail HTML to a file
    def output_d = new File("${params.outdir}/pipeline_info/")
    if (!output_d.exists()) {
        output_d.mkdirs()
    }
    def output_hf = new File(output_d, "${runName}_pipeline_report.html")
    output_hf.withWriter { w -> w << email_html }
    def output_tf = new File(output_d, "${runName}_pipeline_report.txt")
    output_tf.withWriter { w -> w << email_txt }

    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_red = params.monochrome_logs ? '' : "\033[0;31m";
    c_reset = params.monochrome_logs ? '' : "\033[0m";

    if (workflow.stats.ignoredCount > 0 && workflow.success) {
        log.info "-${c_purple}Warning, pipeline completed, but with errored process(es) ${c_reset}-"
        log.info "-${c_red}Number of ignored errored process(es) : ${workflow.stats.ignoredCount} ${c_reset}-"
        log.info "-${c_green}Number of successfully ran process(es) : ${workflow.stats.succeedCount} ${c_reset}-"
    }

    if (workflow.success) {
        log.info "-${c_purple}[jsimonas/solo-in-drops]${c_green} Pipeline completed successfully${c_reset}-"
    } else {
        checkHostname()
        log.info "-${c_purple}[jsimonas/solo-in-drops]${c_red} Pipeline completed with errors${c_reset}-"
    }

}

def nfcoreHeader() {
    // Log colors ANSI codes
    c_black = params.monochrome_logs ? '' : "\033[0;30m";
    c_blue = params.monochrome_logs ? '' : "\033[0;34m";
    c_cyan = params.monochrome_logs ? '' : "\033[0;36m";
    c_dim = params.monochrome_logs ? '' : "\033[2m";
    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_reset = params.monochrome_logs ? '' : "\033[0m";
    c_white = params.monochrome_logs ? '' : "\033[0;37m";
    c_yellow = params.monochrome_logs ? '' : "\033[0;33m";

    return """    -${c_dim}--------------------------------------------------${c_reset}-
                                            ${c_green},--.${c_black}/${c_green},-.${c_reset}
    ${c_blue}        ___     __   __   __   ___     ${c_green}/,-._.--~\'${c_reset}
    ${c_blue}  |\\ | |__  __ /  ` /  \\ |__) |__         ${c_yellow}}  {${c_reset}
    ${c_blue}  | \\| |       \\__, \\__/ |  \\ |___     ${c_green}\\`-._,-`-,${c_reset}
                                            ${c_green}`._,._,\'${c_reset}
    ${c_purple}  jsimonas/solo-in-drops v${workflow.manifest.version}${c_reset}
    -${c_dim}--------------------------------------------------${c_reset}-
    """.stripIndent()
}

def checkHostname() {
    def c_reset = params.monochrome_logs ? '' : "\033[0m"
    def c_white = params.monochrome_logs ? '' : "\033[0;37m"
    def c_red = params.monochrome_logs ? '' : "\033[1;91m"
    def c_yellow_bold = params.monochrome_logs ? '' : "\033[1;93m"
    if (params.hostnames) {
        def hostname = "hostname".execute().text.trim()
        params.hostnames.each { prof, hnames ->
            hnames.each { hname ->
                if (hostname.contains(hname) && !workflow.profile.contains(prof)) {
                    log.error "====================================================\n" +
                            "  ${c_red}WARNING!${c_reset} You are running with `-profile $workflow.profile`\n" +
                            "  but your machine hostname is ${c_white}'$hostname'${c_reset}\n" +
                            "  ${c_yellow_bold}It's highly recommended that you use `-profile $prof${c_reset}`\n" +
                            "============================================================"
                }
            }
        }
    }
}