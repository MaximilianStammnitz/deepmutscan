//
// Subworkflow with functionality specific to the nf-core/deepmutscan pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFSCHEMA_PLUGIN     } from '../../nf-core/utils_nfschema_plugin'
include { paramsSummaryMap          } from 'plugin/nf-schema'
include { samplesheetToList         } from 'plugin/nf-schema'
include { completionEmail           } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_INITIALISATION {

    take:
    version           // boolean: Display version and exit
    validate_params   // boolean: Boolean whether to validate parameters against the schema at runtime
    monochrome_logs   // boolean: Do not use coloured log outputs
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir            //  string: The output directory where the results will be saved
    input             //  string: Path to input samplesheet

    main:

    ch_versions = Channel.empty()

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE (
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1
    )

    //
    // Validate parameters and generate parameter summary to stdout
    //

    def before_text = ""
    def after_text = ""
    before_text = """
-\033[2m----------------------------------------------------\033[0m-
                                        \033[0;32m,--.\033[0;30m/\033[0;32m,-.\033[0m
\033[0;34m        ___     __   __   __   ___     \033[0;32m/,-._.--~\'\033[0m
\033[0;34m  |\\ | |__  __ /  ` /  \\ |__) |__         \033[0;33m}  {\033[0m
\033[0;34m  | \\| |       \\__, \\__/ |  \\ |___     \033[0;32m\\`-._,-`-,\033[0m
                                        \033[0;32m`._,._,\'\033[0m
\033[0;35m  nf-core/deepmutscan ${workflow.manifest.version}\033[0m
-\033[2m----------------------------------------------------\033[0m-
"""
    after_text = """${workflow.manifest.doi ? "\n* The pipeline\n" : ""}${workflow.manifest.doi.tokenize(",").collect { doi -> "    https://doi.org/${doi.trim().replace('https://doi.org/','')}"}.join("\n")}${workflow.manifest.doi ? "\n" : ""}
* The nf-core framework
    https://doi.org/10.1038/s41587-020-0439-x

* Software dependencies
    https://github.com/nf-core/deepmutscan/blob/master/CITATIONS.md
"""
    if (monochrome_logs) {
        before_text = before_text.replaceAll(/\033\[[0-9;]*m/, '')
    }

    command = "nextflow run ${workflow.manifest.name} -profile <docker/singularity/.../institute> --input samplesheet.csv --outdir <OUTDIR>"

    UTILS_NFSCHEMA_PLUGIN (
        workflow,
        validate_params,
        params.help,
        params.help_full,
        params.show_hidden,
        null,
        before_text,
        after_text,
        command,
        false
    )

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE (
        nextflow_cli_args
    )

    //
    // Custom validation for pipeline parameters
    //
    validateInputParameters()

    //
    // Create channel from input file provided through params.input
    //
    // samplesheetToList() reads the samplesheet with assets/schema_input.json: every row is
    // validated (FASTQ extensions, both reads present, allowed types, no duplicate
    // sample/type/replicate) and file paths are returned as file objects.
    //
    def samplesheet_rows = samplesheetToList(input, "${projectDir}/assets/schema_input.json")
        .collect { meta, file1, file2 ->
            def new_meta = [
                id        : "${meta.sample}_${meta.type}_${meta.replicate}_pe".toString(),
                sample    : meta.sample,
                type      : meta.type,
                replicate : meta.replicate as int
            ]
            [ new_meta, [ file1, file2 ] ]
        }

    validateSamplesheet(samplesheet_rows.collect { row -> row[0] })

    ch_samplesheet = channel.fromList(samplesheet_rows)

    // Emit the samplesheet channel and an empty version channel for use in the workflow
    emit:
    samplesheet = ch_samplesheet
    versions    = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_COMPLETION {

    take:
    email           //  string: email address
    email_on_fail   //  string: email address sent on pipeline failure
    plaintext_email // boolean: Send plain-text email instead of HTML
    outdir          //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output
    multiqc_report  //  string: Path to MultiQC report

    main:
    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def multiqc_reports = multiqc_report.toList()

    //
    // Completion email and summary
    //
    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(
                summary_params,
                email,
                email_on_fail,
                plaintext_email,
                outdir,
                monochrome_logs,
                multiqc_reports.getVal(),
            )
        }

        completionSummary(monochrome_logs)

    }

    workflow.onError {
        log.error "Pipeline failed. Please refer to troubleshooting docs for common issues: https://nf-co.re/docs/running/troubleshooting"
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
//
// Check and validate pipeline parameters
//
def validateInputParameters() {
    genomeExistsError()
    readingFrameError()
    customCodonLibraryError()
}

//
// Check that --reading_frame covers whole codons and lies within the reference FASTA
//
def readingFrameError() {
    if (!params.reading_frame || !params.fasta) {
        return  // both are required parameters; the schema validation reports them if missing
    }
    def (start, stop) = params.reading_frame.tokenize('-').collect { coord -> coord as long }
    if (start < 1 || stop < start) {
        error("[--reading_frame ${params.reading_frame}] The reading frame must be given as 'start-stop' with 1 <= start <= stop (1-based, inclusive).")
    }
    def length = stop - start + 1
    if (length % 3 != 0) {
        error("[--reading_frame ${params.reading_frame}] The reading frame spans ${length} nucleotides, which is not a multiple of 3. 'start' must be the first nucleotide of the first codon and 'stop' the last nucleotide of the last codon.")
    }
    def fasta_length = fastaSequenceLength(params.fasta)
    if (stop > fasta_length) {
        error("[--reading_frame ${params.reading_frame}] The reading frame ends at nucleotide ${stop}, but the reference sequence in '--fasta' is only ${fasta_length} nucleotides long.")
    }
}

//
// Length of the (first) sequence in a FASTA file, plain or gzipped
//
def fastaSequenceLength(fasta_path) {
    def fasta  = file(fasta_path, checkIfExists: true)
    def stream = fasta.name.endsWith('.gz') ? new java.util.zip.GZIPInputStream(fasta.newInputStream()) : fasta.newInputStream()
    def length = 0L
    def n_records = 0
    stream.withReader { reader ->
        reader.eachLine { line ->
            if (line.startsWith('>')) {
                n_records += 1
            } else if (n_records == 1) {
                length += line.trim().length()
            }
        }
    }
    return length
}

//
// '--mutagenesis_type custom' needs a custom codon library
//
def customCodonLibraryError() {
    if (params.mutagenesis_type == 'custom' && (!params.custom_codon_library || params.custom_codon_library == '/NULL')) {
        error("[--mutagenesis_type custom] Please provide the custom codon library with '--custom_codon_library <file.csv>'.")
    }
}

//
// Check the samplesheet as a whole (row-level checks are done by assets/schema_input.json)
//
def validateSamplesheet(metas) {
    // One biological sample per run
    def samples = metas.collect { meta -> meta.sample }.unique()
    if (samples.size() > 1) {
        error("[samplesheet] nf-core/deepmutscan v1.0.0 supports one biological sample per run (found: ${samples.join(', ')}). Please split the samplesheet and run the pipeline once per sample.")
    }

    // Fitness estimation pairs input and output libraries by replicate number
    if (params.fitness) {
        def input_reps  = metas.findAll { meta -> meta.type == 'input'  }.collect { meta -> meta.replicate }.sort()
        def output_reps = metas.findAll { meta -> meta.type == 'output' }.collect { meta -> meta.replicate }.sort()
        if (!input_reps || !output_reps) {
            error("[--fitness] Fitness estimation needs at least one 'input' and one 'output' library in the samplesheet (found ${input_reps.size()} input and ${output_reps.size()} output).")
        }
        if (input_reps != output_reps) {
            error("[--fitness] The 'input' and 'output' libraries must have the same replicate numbers, as they are paired by replicate (input: ${input_reps.join(', ')}; output: ${output_reps.join(', ')}).")
        }
    }

    // '--error_correction wildtype' needs exactly one wildtype library
    if (params.error_correction == 'wildtype') {
        def n_wildtype = metas.count { meta -> meta.type == 'wildtype' }
        if (n_wildtype == 0) {
            error(
                "[--error_correction wildtype] No 'wildtype' library found in the samplesheet.\n" +
                "  Add a samplesheet row with type 'wildtype' and sample '${samples[0]}' (deep wildtype-only sequencing),\n" +
                "  or use the internal correction instead (default: --error_correction false_doubles),\n" +
                "  or disable error correction with --error_correction none."
            )
        }
        if (n_wildtype > 1) {
            error("[--error_correction wildtype] Found ${n_wildtype} 'wildtype' libraries in the samplesheet; exactly one is supported. Please keep one wildtype library (e.g. the deepest-sequenced one).")
        }
    }
}

//
// Get attribute from genome config file e.g. fasta
//
def getGenomeAttribute(attribute) {
    if (params.genomes && params.genome && params.genomes.containsKey(params.genome)) {
        if (params.genomes[ params.genome ].containsKey(attribute)) {
            return params.genomes[ params.genome ][ attribute ]
        }
    }
    return null
}

//
// Exit pipeline if incorrect --genome key provided
//
def genomeExistsError() {
    if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
        def error_string = "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
            "  Genome '${params.genome}' not found in any config files provided to the pipeline.\n" +
            "  Currently, the available genome keys are:\n" +
            "  ${params.genomes.keySet().join(", ")}\n" +
            "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        error(error_string)
    }
}
//
// Generate methods description for MultiQC
//
def toolCitationText() {
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "Tool (Foo et al. 2023)" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def citation_text = [
            "Tools used in the workflow included:",
            "FastQC (Andrews 2010),",
            "MultiQC (Ewels et al. 2016)",
            "."
        ].join(' ').trim()

    return citation_text
}

def toolBibliographyText() {
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "<li>Author (2023) Pub name, Journal, DOI</li>" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def reference_text = [
            "<li>Andrews S, (2010) FastQC, URL: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).</li>",
            "<li>Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics , 32(19), 3047–3048. doi: /10.1093/bioinformatics/btw354</li>"
        ].join(' ').trim()

    return reference_text
}

def methodsDescriptionText(mqc_methods_yaml) {
    // Convert  to a named map so can be used as with familiar NXF ${workflow} variable syntax in the MultiQC YML file
    def meta = [:]
    meta.workflow = workflow.toMap()
    meta["manifest_map"] = workflow.manifest.toMap()

    // Pipeline DOI
    if (meta.manifest_map.doi) {
        // Using a loop to handle multiple DOIs
        // Removing `https://doi.org/` to handle pipelines using DOIs vs DOI resolvers
        // Removing ` ` since the manifest.doi is a string and not a proper list
        def temp_doi_ref = ""
        def manifest_doi = meta.manifest_map.doi.tokenize(",")
        manifest_doi.each { doi_ref ->
            temp_doi_ref += "(doi: <a href=\'https://doi.org/${doi_ref.replace("https://doi.org/", "").replace(" ", "")}\'>${doi_ref.replace("https://doi.org/", "").replace(" ", "")}</a>), "
        }
        meta["doi_text"] = temp_doi_ref.substring(0, temp_doi_ref.length() - 2)
    } else meta["doi_text"] = ""
    meta["nodoi_text"] = meta.manifest_map.doi ? "" : "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

    // Tool references
    meta["tool_citations"] = ""
    meta["tool_bibliography"] = ""

    meta["tool_citations"] = toolCitationText().replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
    meta["tool_bibliography"] = toolBibliographyText()


    def methods_text = mqc_methods_yaml.text

    def engine =  new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}
