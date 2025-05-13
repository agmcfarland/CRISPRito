// ----------------------------------------
// PARAMETERS
// ----------------------------------------
params.sample_sheet        = "data/input_sample_sheet.csv"
params.genome              = "data/hg38.fa.gz"
params.features            = "data/features.csv"
params.genes               = "data/gene_names.csv"
params.output_dir          = "CRISPRito_output"
params.overwrite_output_dir = true
params.flank_size          = 30
params.sgRNA               = ""
params.pam_alignment       = "-GG"
// params.scripts_dir         = "/data/CRISPRito/CRISPRito"

// ----------------------------------------
// CHANNEL DEFINITIONS
// ----------------------------------------

sample_sheet_ch = Channel.fromPath(params.sample_sheet)
genome_ch       = Channel.fromPath(params.genome)
features_ch     = Channel.fromPath(params.features)
genes_ch        = Channel.fromPath(params.genes)
overwrite_ch    = Channel.value(params.overwrite_output_dir ? "--overwrite_output_dir" : "")

// ----------------------------------------
// WORKFLOW
// ----------------------------------------

workflow {
	setup_out = setup_run(
		sample_sheet_ch,
		genome_ch,
		features_ch,
		genes_ch,
		overwrite_ch
	)

	// Flatten the output channel so process_groups gets one file at a time
	setup_out.out.flatten().set { group_sheet_ch }

	process_groups(group_sheet_ch, genome_ch, features_ch, genes_ch)
}
// ----------------------------------------
// Step 1: Setup Run
// ----------------------------------------

process setup_run {

	publishDir "${params.output_dir}", mode: 'copy'

	input:
	path sample_sheet
	path genome
	path features
	path genes
	val overwrite_flag

	output:
	  path("*_group_samplesheet.csv"), emit: out

	conda 'CRISPRito'

	script:
	"""
	setup_run \\
		--sample_sheet_path ${sample_sheet} \\
		--output_dir ${params.output_dir} \\
		--genome_path ${genome} \\
		--feature_path ${features} \\
		--gene_names_path ${genes} \\
		${overwrite_flag}

	cp ${params.output_dir}/*_group_samplesheet.csv .
	"""
}

// ----------------------------------------
// Step 2: Process Groups
// ----------------------------------------

process process_groups {

	publishDir "${params.output_dir}", mode: 'copy'

	input:
	path group_sheet
	path genome
	path features
	path genes

	output:
	path("*_group_*.csv")

	conda 'CRISPRito'

	script:
	"""
	process_group \\
		${group_sheet} \\
		--output_dir ${params.output_dir} \\
		--genome_path ${genome} \\
		--feature_path ${features} \\
		--gene_names_path ${genes} \\
		--flank_size ${params.flank_size} \\
		--sgRNA ${params.sgRNA} \\
		--PAM_alignment ${params.pam_alignment}
	"""
}

