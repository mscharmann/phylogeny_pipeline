configfile: "config.yaml"

genomefile = config["genomefile"]
samples_units_fqs_map = config["samples_units_fqs_map"]
gff = config["gff"]

import pandas as pd

samples_units_fqs = pd.read_table(samples_units_fqs_map, dtype=str).set_index(
	["sample", "unit", "fq1", "fq2"], drop=False)

SAMPLES = list( set(samples_units_fqs["sample"]) )


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


def read_wrapped_or_unwrapped_fasta (infile):
	outlines = []
	with open(infile, "r") as INFILE:
		first_id = INFILE.readline().strip("\n").strip(">")
		outlines.append(first_id)
		seq = ""
		for line in INFILE:
			line = line.strip("\n")
			if line.startswith(">"):
				outlines.append(seq)
				outlines.append(line.strip(">").strip("\n"))
				seq = ""
			else:
				if len(line) > 0:
					seq += line.strip("\n").upper()
		# append last seq
		outlines.append(seq)
	i=0
	j=1
	out_dict = {}
	for x in range(int(len(outlines)/2)):
		out_dict[outlines[i]] = outlines[j]
		i += 2
		j += 2
	return out_dict



rule all:
	input:
		expand("mapped_reads/{sample}.sorted.bam.bai", sample=SAMPLES),
		expand("results_raw/features.{sample}.fa", sample=SAMPLES),
		"results_processed/supermatrix.stats.txt",
		"results_processed/mapping_statistics_report.txt"
	run:
		import os

		os.system( "rm -r varcall*" )
		os.system( "rm -r mapped_reads_per_unit" )


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
		sleep 20
		tabix {output}
		"""



rule VCF_filter_variants:
	input:
		gzvcf="varcall_chunk_VCFs/{i}.bed.vcf.gz"
	output:
		temp( "varcall_chunk_VCFs_filtered/{i}.bed.vcf.gz" )
	params:
		QUAL=config["VCF_QUAL"],
		MIN_DEPTH=config["MIN_DEPTH"]
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

		# variants: || "vcftools --min-meanDP" is a bad idea here, because it is a site-filter, not a genotype-filter: some low-coverage samples in the set may cause mean depth across samples to be low, thus loosing sites that are good and credible in the other samples!
		#vcftools --gzvcf ../../{input.gzvcf} --mac 1 --minQ {params.QUAL} --min-meanDP {params.MIN_DEPTH} --minDP {params.MIN_DEPTH} --max-meanDP $MAX_DEPTH --recode --stdout | bgzip -c > tmp.1
		vcftools --gzvcf ../../{input.gzvcf} --mac 1 --minQ {params.QUAL} --minDP {params.MIN_DEPTH} --max-meanDP $MAX_DEPTH --recode --stdout | bgzip -c > tmp.1
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
		temp("results_raw/low_cov_regions.{sample}.bed")
	params:
		MIN_DEPTH=config["MIN_DEPTH"]
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
		bed=temp("results_raw/high_cov_regions.{sample}.bed"),
		d=temp("results_raw/depth.{sample}.txt"),
		dc="results_raw/depth_cutoff.{sample}.txt"
	shell:
		"""
		samtools depth {input.bam} | awk 'NR % 10000 == 0 {{print $3}}' > {output.d}

		cat {output.d} | sort -n | awk '{{all[NR] = $0}} END{{print all[int(NR*0.96 - 0.5)]}}' > {output.dc}

		cut=$(cat {output.dc}	)

		# set an absolute minimum for the max depth cutoff: with really low-coverage data, such as "genome skimming", the 96% percentile is often around 3 reads.
		# Thus increase the maximum depth cutoff to something like 15.

		cut=$(cat {output.dc} )
		if [[ $cut -lt 15 ]]; then
			cut=15
			echo $cut > {output.dc}
		fi

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


rule apply_variants_and_gaps_to_genome_and_extract_features:
	input:
		bad="results_raw/cov_excluded_regions.{sample}.bed",
		vcf="results_raw/variants.post_filter.vcf.gz",
		fa=genomefile,
		gff=gff
	output:
		temp( "results_raw/features.{sample}.fa" )
	shell:
		"""
		sname=$( echo {output} | sed 's/results_raw\/features\.//g' | sed 's/\.fa//g' )
		bcftools consensus -f {input.fa} -m {input.bad} --haplotype I --missing N -s $sname {input.vcf} > $sname.GENOME_applied_variants_and_gaps.fa
		sleep 10
		gffread -g $sname.GENOME_applied_variants_and_gaps.fa -x {output} {input.gff}
		sleep 2
		rm $sname.GENOME_applied_variants_and_gaps.fa $sname.GENOME_applied_variants_and_gaps.fa.fai
		"""


rule bin_per_feature:  # this rule takes fastas with features per sample and re-arranges to one fasta per feature
	input:
		expand("results_raw/features.{sample}.fa", sample=SAMPLES)
	output:
		temp( "feature_fastas_success" )
	run:
		# python code
		import sys, os

		os.system("rm -r feature_fastas")

		allfiles = os.listdir("results_raw")
		infiles = sorted([x for x in allfiles if x.startswith("features.")])

		print("parsing files to get all features")

		all_features = set()
		for f in infiles:
			with open("results_raw/" + f, "r") as I:
				for line in I:
					if line.startswith(">"):
						t = line.strip(">").strip()
						all_features.add(t)

		all_features = list(sorted(all_features))
		print("features:")
		print (len(all_features))

		print("second pass, now to get sequences and assembling outputs")

		indata_dict = {}
		for f in infiles:
			inseqs = read_wrapped_or_unwrapped_fasta ("results_raw/" + f)
			indata_dict[f] = inseqs

		os.mkdir("feature_fastas")

		# split list into chunks of size n
		n = 1000
		x = [all_features[i:i + n] for i in range(0, len(all_features), n)]

		cnt = 0
		for featurechunk in x:
			cnt += 1
			print(cnt)
			for feature in featurechunk:
				outlines = []
				for f in infiles:
					outlines.append(">" + ".".join(f.split(".")[1:-1]))
					outlines.append( indata_dict[f][feature] )
				with open("feature_fastas/" + feature + ".fa", "a") as O:
					O.write("\n".join(outlines)+"\n")
		os.system("touch feature_fastas_success")



rule filter_column_missingness:
	input:
		"feature_fastas_success"
	output:
		temp( "feature_fastas_filter_success" )
	params:
		COLUMN_MIN_PRESENCE=config["COLUMN_MIN_PRESENCE"],
		MIN_BASES_PER_SAMPLE=config["MIN_BASES_PER_SAMPLE"]
	run:
		# filter column missingness and retain only COMPLETE codons (given input codons!)
		# exclude columns if they are present in < nsam * COLUMN_MIN_PRESENCE
		# exclude samples that have < MIN_BASES_PER_SAMPLE
		# replace N with - character

		column_min_pres_thresh = float( params.COLUMN_MIN_PRESENCE )
		min_bases_per_sample = int( params.MIN_BASES_PER_SAMPLE )

		import os
		os.system("rm -r feature_fastas_missingness_filtered")
		os.system("mkdir feature_fastas_missingness_filtered")

		allfiles = os.listdir("feature_fastas/")
		infiles = sorted([x for x in allfiles if x.endswith(".fa")])

		for f in infiles:
			indict = read_wrapped_or_unwrapped_fasta ("feature_fastas/" + f)
			sample_order = list( indict.keys() )
			nsam = len(sample_order)
			output_columns = []
			out_dict = {k:"" for k in sample_order}
			# problem: sometimes gene predictions contain errors, so the full length of a coding sequence is not a multiple of three. In that case the trailing 1 or 2 sites shall be omitted.
			iter_lim = len( list(indict.values())[0] )
			while iter_lim % 3 != 0:
				iter_lim -= 1
			idx=0
			while idx < iter_lim:
				codon_column = [ indict[x][idx:(idx+3)] for x in sample_order ]
				cleaned_codons = [ x if not "N" in x else "---" for x in codon_column]
				present_codons = len( [ x for x in cleaned_codons if not "-" in x ] )
				if present_codons >= column_min_pres_thresh*nsam:
					for sidx,s in enumerate(sample_order):
						out_dict[s] += cleaned_codons[sidx]
				idx += 3
			final_out_dict = {}
			for k,v in out_dict.items():
				if len( [ x for x in v if x != "-" ] ) >= min_bases_per_sample:
					final_out_dict[k] = v
			outlines = [">" + k + "\n" + v for k,v in final_out_dict.items()]
			with open("feature_fastas_missingness_filtered/" + f, "w") as O:
				O.write("\n".join(outlines) + "\n")
		os.system("rm -r feature_fastas")
		os.system("touch feature_fastas_filter_success")




rule assemble_supermatrix:
	input:
		"feature_fastas_filter_success"
	output:
		"results_processed/supermatrix.stats.txt"
	params:
		MIN_LEN_LOCUS = config["MIN_LEN_LOCUS"],
		MIN_SAMPLE_FRACTION = config["MIN_SAMPLE_FRACTION"]
	run:
		# reads fasta alignments per-locus from a directory and assemble a supermatrix
		# excludes alignments that are < MIN_LEN_LOCUS
		# excludes alignments that have < MIN_SAMPLE_FRACTION * n_samples

		import sys, os

		os.system("rm -r results_processed")
		os.system("mkdir results_processed")

		allfiles = os.listdir( "feature_fastas_missingness_filtered" )
		infiles = sorted([x for x in allfiles if x.endswith(".fa")])

		min_len_locus = float( params.MIN_LEN_LOCUS )
		min_sample_fraction = float( params.MIN_SAMPLE_FRACTION )

		print("number of input files:")
		print(len(infiles))
		print("parsing files to get all taxa")

		all_samples = set()
		for f in infiles:
			with open("feature_fastas_missingness_filtered/" + f, "r") as I:
				for line in I:
					if line.startswith(">"):
						t = line.strip(">").strip()
						all_samples.add(t)

		all_samples = list(sorted(all_samples))
		print("n_samples:")
		print (len(all_samples))
		print("second pass, now to get sequences and assembling outputs")

		locus_count = 0
		out_seq_dict = {t:[] for t in all_samples}
		modelfile_lines = []
		left_idx = 1
		right_idx = 0

		for f in infiles:
			nseq = 0
			samples = []
			seqlen = 0
			indict = read_wrapped_or_unwrapped_fasta ("feature_fastas_missingness_filtered/" + f)
			seqlen = len(  list(indict.values())[0] )
			samples = list(indict.keys())
			nseq = len(samples)
			if nseq >= min_sample_fraction*len(all_samples): # sufficient number of samples have this locus
				if seqlen >= min_len_locus:	# locus has sufficient length
					locus_count += 1
					# prep supermatrix; filling in gaps for samples who do not have any data at this locus
					for t in all_samples:
						try:
							out_seq_dict[t].append( indict[t] )
						except KeyError:
							out_seq_dict[t].append( "-"*seqlen )
					# prep modelfile
					right_idx += seqlen
					mline = "DNA" + "," + ".".join(f.split(".")[:-1]) + "=" + str( left_idx ) + "-" + str( right_idx )
					left_idx += seqlen
					modelfile_lines.append( mline )


		with open("results_processed/supermatrix.model.txt","w") as outfile2:
				outfile2.write( "\n".join( modelfile_lines ) + "\n")


		outlines_fasta = [">" + k + "\n" + "".join(v) for k,v in out_seq_dict.items()]
		with open("results_processed/supermatrix.fasta","w") as outfile1:
				outfile1.write( "\n".join( outlines_fasta ) + "\n")


		with open("results_processed/supermatrix.stats.txt", "w") as O:
			O.write("\t".join( [ "sample" , "n_loci" , "total_length" , "proportion_gap" ] ) + "\n")
			for k,v in out_seq_dict.items():
				sline = []
				sline.append(k)
				sline.append( str(len(v)) )
				all = "".join(v)
				sline.append( str( len(all) ) )
				nongap = "".join(v).replace("-","")
				sline.append( str( round( 1.0 - float(len(nongap))/float(len(all)),3  ) ) )
				O.write("\t".join( sline ) + "\n")
			O.write("\n")
			O.write("after filtering exported " + str(locus_count) + " out of " + str(len(infiles)) + " loci" + "\n")
		os.system("rm -r feature_fastas_missingness_filtered")


rule count_reads_input_and_mapping:
	input:
		indexfile="mapped_reads/{sample}.sorted.bam.bai",
		bamfile="mapped_reads/{sample}.sorted.bam",
		reads=get_fastq_sample_ALL
	output:
		temp( "results_raw/{sample}.raw_and_mapped_reads_report.txt" )
	threads: 3
	run:
		cmd = "samtools flagstat " + input.bamfile + " > " + input.bamfile + ".flagstat.txt"
		os.system( cmd )
		with open(input.bamfile +".flagstat.txt", "r") as I:
			total_mapped = int(I.readline().split()[0])
		os.system( "rm " +  input.bamfile +".flagstat.txt")
		# reconstruct the sample name...
		samplename = input.bamfile.split("/")[1].split(".sorted.bam")[0]
		# now count all the reads
		sum_of_read_fastq_lines = 0
		for f in list(input.reads):
			os.system( "pigz -p3 -dc {0} | wc -l > tmp.{1}.readcount".format(f,samplename) )
			with open( "tmp.{}.readcount".format(samplename), "r" ) as an_infile:
				sum_of_read_fastq_lines += int(an_infile.readline().strip())
			os.system( "rm tmp.{}.readcount".format(samplename) )
		n_reads = float(sum_of_read_fastq_lines)/4.0
		mapping_rate = float(total_mapped)/n_reads
		with open(str(output), "w") as O:
			O.write(samplename + "\t" + str(n_reads) + "\t" + str(total_mapped) + "\t" + str(round(mapping_rate,3)) + "\n")


rule collect_mapping_report:
	input:
		expand("results_raw/{sample}.raw_and_mapped_reads_report.txt", sample=SAMPLES)
	output:
		"results_processed/mapping_statistics_report.txt"
	shell:
		"""
		echo -e "counts are single reads, not pairs of reads" > {output}
		echo -e "sample\traw_reads\tmapped_reads\tmapping_rate" >> {output}
		cat results_raw/*.raw_and_mapped_reads_report.txt >> {output}
		"""
