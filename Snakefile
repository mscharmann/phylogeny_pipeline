configfile: "config.yaml"

genomefile = config["genomefile"]
samples_units_fqs_map = config["samples_units_fqs_map"]
windowfile = config["windows_for_gene_trees"]

import pandas as pd

samples_units_fqs = pd.read_table(samples_units_fqs_map, dtype=str).set_index(
	["sample", "unit", "fq1", "fq2"], drop=False)

SAMPLES = list( set(samples_units_fqs["sample"]) ) 
print(samples_units_fqs)


def get_sample_bams(wildcards):
	"""Get all aligned reads of given sample."""
	return expand(
		"mapped_reads_per_unit/{sample}-{unit}.sorted.bam",
		sample=wildcards.sample,
		unit=samples_units_fqs.loc[wildcards.sample].unit,
	)


def get_fastq_sample_unit(wildcards):
	"""Get fastq files of given sample-unit."""
	fastqs = samples_units_fqs.loc[(wildcards.sample, wildcards.unit), ["fq1", "fq2"]]
	if fastqs.fq2.isnull().values.any():
		return [ fastqs.fq1.item() ]
	return [ fastqs.fq1.item(), fastqs.fq2.item() ]


def get_fastq_sample_ALL(wildcards):
	"""Get list of fastq files of given sample including ALL units."""
	fastqs_pd = samples_units_fqs.loc[(wildcards.sample), ["fq1", "fq2"]]
	fastqs = set( fastqs_pd.fq1.tolist() + fastqs_pd.fq2.tolist() )
	fastqs_clean = list( {x for x in fastqs if pd.notna(x)} )
	return fastqs_clean



rule all:
	input:
		expand("mapped_reads/{sample}.sorted.bam.bai", sample=SAMPLES),
		expand("results_raw/cov_excluded_regions.{sample}.bed", sample=SAMPLES),
		"results_raw/variants.post_filter.vcf.gz",
		"window_fastas"
	shell:
		"""
		rm -r varcall*
		rm -r mapped_reads_per_unit
		"""



rule bwa_idx:
	input:
	   genomefile
	output:
		"{genomefile}.bwt"
	shell:
		"""
		if [[ ! $( grep ">" {input} ) =~ "|" ]]; then
			bwa index {input}
		else
			echo "refusing to run, fasta headers contain pipe '|' character, dying"
		fi
		"""

rule bwa_map:
	input:
		fa=genomefile,
		gidx=genomefile + ".bwt",
		reads=get_fastq_sample_unit
	output:
		temp("mapped_reads_per_unit/{sample}-{unit}.bam")
	threads: 12
	run:
		if len(input.reads) == 2: # paired-end!
			shell("""
				# filtering alignments to be primary (i.e. each read only once; if multiple locations equally possible than a random one is chosen): 
				# -F 256 == -F 0x0100 == NOT not primary alignment
				# filtering alignments to be NOT supplementary (supplementary: sections of read map to discontinuous coordinates, e.g. across an inversion breakpoint..):
				# -F 2048 == -F 0x800 == NOT supplementary alignment
				# sum of the bit flags: 2304 => filters against BOTH non-primary and supplementary alignments; verified with samtools flagstat
				# filtering alignments to be "properly paired": -f 2
				# filtering against multi-mapping alignments (which have MAPQ=0): -q 1
				bwa mem -t {threads} -a {input.fa} {input.reads[0]} {input.reads[1]} -R "@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}\\tPL:Illumina" | samtools view -F 2304 -f 2 -q 1 -b -@ 2 - > {output} 
				""")
		else: # single-end
			shell("""
				# filtering alignments to be primary (i.e. each read only once; if multiple locations equally possible than a random one is chosen): 
				# -F 256 == -F 0x0100 == NOT not primary alignment
				# filtering alignments to be NOT supplementary (supplementary: sections of read map to discontinuous coordinates, e.g. across an inversion breakpoint..):
				# -F 2048 == -F 0x800 == NOT supplementary alignment
				# -F 4 read unmapped (0x4)
				# sum of the bit flags: 2308 => filters against non-primary and supplementary alignments and unmapped
				# filtering against multi-mapping alignments (which have MAPQ=0): -q 1
				bwa mem -t {threads} -a {input.fa} {input.reads[0]} -R "@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}\\tPL:Illumina" | samtools view -F 2308 -q 1 -b -@ 2 - > {output} 
				""")


rule samtools_sort:
	input:
		"mapped_reads_per_unit/{sample}-{unit}.bam"
	output:
		temp( "mapped_reads_per_unit/{sample}-{unit}.sorted.bam" )
	shell:
		"""
		samtools sort -T mapped_reads_per_unit/{wildcards.sample}.{wildcards.unit} -O bam {input} > {output}
		"""
		
rule merge_bams_per_sample:
	input:
		bams=get_sample_bams
	output:
		"mapped_reads/{sample}.sorted.bam"
	threads: 4
	shell:
		"""
		samtools merge --threads {threads} {output} {input.bams}
		"""


rule samtools_index:
	input:
		"mapped_reads/{sample}.sorted.bam"
	output:
		"mapped_reads/{sample}.sorted.bam.bai"
	shell:
		"samtools index {input}"




# A checkpoint that shall trigger re-evaluation of the DAG
# https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#data-dependent-conditional-execution
# we need this because the number of output files (i.e. chunks / region files for variant calling) is not known before execution!
checkpoint split_ref_for_varcall:
	input:
		ref=genomefile,
		indexes=expand("mapped_reads/{sample}.sorted.bam.bai", sample=SAMPLES),
		samples=expand("mapped_reads/{sample}.sorted.bam", sample=SAMPLES)
	output:
		directory("varcall_chunks")
	params:
		chunksize_Mb = config["varcall_chunksize"]
	shell:
		"""
		samtools faidx {input.ref}
		
		echo {input.samples} | sed 's/ /\\n/g' > bamlistforsplit
		# wget https://gist.githubusercontent.com/travc/0c53df2c8eca81c3ebc36616869930ec/raw/eff3032ca7c955ca33bffd8758092e4006949c75/split_ref_by_bai_datasize.py
		python scripts/split_ref_by_bai_datasize.py -r {input.ref}.fai -L bamlistforsplit --target-data-size {params.chunksize_Mb} > varcall_regionsALL.bed
		rm bamlistforsplit
		split --lines=1 varcall_regionsALL.bed varcall_regions_ --numeric-suffixes --suffix-length=6 --additional-suffix=.bed
		mkdir -p {output}
		for i in varcall_regions_* ; do mv $i {output} ; done		
		"""

rule call_variants:
	input:
		ref=genomefile,
		regions="varcall_chunks/{i}.bed"
	output:
		temp( "varcall_chunk_VCFs/{i}.bed.vcf.gz" )
	shell:
		"""
		
		# --max-depth 250: use at most 250 reads per input BAM file, apparently these are sampled RANDOMLY!?
		# --min-MQ 15: minimum mapping quality of an alignment, otherwise skip
		# --no-BAQ : do NOT re-calculate mapping quality (which involves re-aligning). Instead, will use MAPQ as stored in BAM file.
		# --min-BQ INT        skip bases with baseQ/BAQ smaller than INT [13]
		# --variants-only
		
		bcftools mpileup -Ou -f {input.ref} -R {input.regions} --bam-list <( ls mapped_reads/*.bam ) --max-depth 250 --min-MQ 20 --min-BQ 15 --no-BAQ -a DP | bcftools call -m --variants-only --skip-variants indels -Ov | bgzip -c > {output}
		
		"""


def merge_vcfs_input(wildcards):
	checkpoint_output = checkpoints.split_ref_for_varcall.get(**wildcards).output[0]
	# print(checkpoint_output)
	vcfs_before_filter = sorted( expand("varcall_chunk_VCFs/{i}.bed.vcf.gz", i=glob_wildcards(os.path.join(checkpoint_output, "{i}.bed")).i) )
	# print (thef)
	vcfs_after_filter = sorted( expand("varcall_chunk_VCFs_filtered/{i}.bed.vcf.gz", i=glob_wildcards(os.path.join(checkpoint_output, "{i}.bed")).i) )
	return vcfs_after_filter

                      

rule merge_filtered_vcfs:
	input:
		region_vcfs=merge_vcfs_input # refers to the function above which evaluates the checkpoint
	output:
		"results_raw/variants.post_filter.vcf.gz"
	threads: 3
	shell:
		"""
		# MUST NOT USE bcftools concat: it cannot resolve POS that are non-monotonically icreasing (which ca happen at the interval boundaries)
		# the code below is only slightly modified from freebayes-parallel script: https://github.com/freebayes/freebayes/blob/master/scripts/freebayes-parallel
		# zcat input.region_vcfs | python $(which vcffirstheader) | vcfstreamsort -w 10000 | vcfuniq | bgzip -c > {output}
		# zcat alone may complain about too many arguments, so better use find -exec :
		find varcall_chunk_VCFs_filtered/*.bed.vcf.gz -type f -exec zcat {{}} \\; | python $(which vcffirstheader) | vcfstreamsort -w 10000 | vcfuniq | bgzip -c > {output}
		tabix {output} 
		"""


	
rule VCF_filter_variants:
	input:
		gzvcf="varcall_chunk_VCFs/{i}.bed.vcf.gz"
	output:
		temp( "varcall_chunk_VCFs_filtered/{i}.bed.vcf.gz" )
	params:
		QUAL=config["VCF_QUAL"],
		MIN_DEPTH=config["VCF_MIN_DEPTH"]
	shell:
		"""
		## https://bcbio.wordpress.com/2013/10/21/updated-comparison-of-variant-detection-methods-ensemble-freebayes-and-minimal-bam-preparation-pipelines/#comment-1469
		## => "We use ... a very vanilla filter with only depth and quality (QUAL < 20, DP < 5)"
		## => impose a minimum of QUAL, and a minimum of DEPTH
		## 	for variants, I want those two, AND a maximum DPETH cutoff 
		## 	for invariants, we can only impose a minimum DEPTH and a maximum DEPTH.
		## use code from dDocent to calcualte a mean depth histogram, then find the 95th percentile: no. too complicated; there are easier ways to do this.
		## also good info:
		## 	https://speciationgenomics.github.io/filtering_vcfs/
		
		# prepare
		wd=DIR_{wildcards.i}
		mkdir -p varcall_chunk_VCFs_filtered/$wd
		cd varcall_chunk_VCFs_filtered/$wd
		
		# set filters
		# MAX_DEPTH is hardcoded here; it is just to reduce the VCF size already at this step. Actual high-coverage filtering is done LATER using regions in .BED format
		MAX_DEPTH=250
			
		# variants:
		vcftools --gzvcf ../../{input.gzvcf} --mac 1 --minQ {params.QUAL} --min-meanDP {params.MIN_DEPTH} --minDP {params.MIN_DEPTH} --max-meanDP $MAX_DEPTH --recode --stdout | bgzip -c > tmp.1
		tabix tmp.1
		
		# filter for missingness: keep any site with at least two samples present
		# first find out how many samples there are:
		nsamples=$( bcftools query -l tmp.1 | wc -l )
		echo $nsamples
		total_minus_two=$((nsamples-2))
		echo $total_minus_two
		vcftools --gzvcf tmp.1 --max-missing-count $total_minus_two --mac 1 --recode --stdout | bcftools view --exclude-uncalled --trim-alt-alleles | bgzip -c > ../../{output}
		
		# cleanup
		cd ../
		rm -r $wd
		
		"""


rule find_low_coverage_regions:
	input:
		bamindex="mapped_reads/{sample}.sorted.bam.bai",
		bam="mapped_reads/{sample}.sorted.bam",
		fa=genomefile
	output:
		"results_raw/low_cov_regions.{sample}.bed"
	params:
		MIN_DEPTH=config["VCF_MIN_DEPTH"]
	shell:
		"""
		genomeCoverageBed -ibam {input.bam} -g {input.fa} -bga -split | awk '{{ if ($4<{params.MIN_DEPTH}) print}}' > {output}
		bedtools merge -i {output} > {output}.tmp
		sleep 2
		mv {output}.tmp {output}

		"""

rule find_high_coverage_regions:
	input:
		bamindex="mapped_reads/{sample}.sorted.bam.bai",
		bam="mapped_reads/{sample}.sorted.bam",
		fa=genomefile
	output:
		bed="results_raw/high_cov_regions.{sample}.bed",
		d=temp("results_raw/depth.{sample}.txt"),
		dc="results_raw/depth_cutoff.{sample}.txt"
	shell:
		"""
		samtools depth {input.bam} | awk 'NR % 10000 == 0 {{print $3}}' > {output.d}

		cat {output.d} | sort -n | awk '{{all[NR] = $0}} END{{print all[int(NR*0.96 - 0.5)]}}' > {output.dc}
		
		cut=$(cat {output.dc}	)
		
		genomeCoverageBed -ibam {input.bam} -g {input.fa} -bga -split | awk -v xfu="$cut" '{{ if ($4>xfu) print}}' > {output.bed}							
		bedtools merge -i {output.bed} > {output.bed}.tmp
		sleep 2
		mv {output.bed}.tmp {output.bed}
		
		"""
		
rule merge_bad_regions:		
	input:
		high="results_raw/high_cov_regions.{sample}.bed",
		low="results_raw/low_cov_regions.{sample}.bed"
	output:
		"results_raw/cov_excluded_regions.{sample}.bed"
	shell:
		"""
		cat <(cut -f 1,2,3 {input.low}) <( cut -f1,2,3 {input.high} ) | bedtools sort > {output}.before_merge

		bedtools merge -i {output}.before_merge > {output}
		rm {output}.before_merge
		"""


rule apply_variants_and_gaps_to_genome:
	input:
		bad="results_raw/cov_excluded_regions.{sample}.bed",
		vcf="results_raw/variants.post_filter.vcf.gz",
		fa=genomefile
	output:
		temp( "results_raw/genome.applied_variants_and_gaps.{sample}.fa" )
	shell:
		"""
		sname=$( echo {output} | sed 's/results_raw\/genome\.applied_variants_and_gaps\.//g' | sed 's/\.fa//g' )
		bcftools consensus -f {input.fa} -m {input.bad} --haplotype I --missing - -s $sname {input.vcf} > {output}
		
		"""


rule get_window_fastas:
	input:
		windows=windowfile,
		g_apps=expand( "results_raw/genome.applied_variants_and_gaps.{sample}.fa", sample=SAMPLES),
		vcf="results_raw/variants.post_filter.vcf.gz"
	output:
		"window_fastas"
	shell:
		"""
		mkdir {output}

		while read w ; do
		wname=$(echo $w | awk '{{print $1"_-_"$2"_-_"$3}}')
		for case in $(bcftools query -l {input.vcf} ); do
			#echo $case
			echo ">"$case > ${{case}}_${{wname}}.temp.fasta 
			echo $w | awk -v var=$case '{{print "samtools faidx results_raw/genome.applied_variants_and_gaps."var".fa "$1":"$2+1"-"$3}}' > fu ; source fu | grep -v ">" | tr -d "\n" | tr 'N' "-" >> ${{case}}_${{wname}}.temp.fasta
			sed -i -e '$a\\' ${{case}}_${{wname}}.temp.fasta
			cat ${{case}}_${{wname}}.temp.fasta >> window_${{wname}}.fasta
			rm ${{case}}_${{wname}}.temp.fasta 
		done
		python2.7 scripts/filter_column_presence.py 0.9 window_${{wname}}.fasta > window_fastas/window_${{wname}}.filtered.fasta
		rm window_${{wname}}.fasta
		done < {input.windows}
		
		rm fu
		"""

		

	
