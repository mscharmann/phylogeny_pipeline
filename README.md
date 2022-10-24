# phylogeny pipeline

inputs:
- reference genome
- GFF file with features named "transcript" (these could be actual genes, or anything else). Each feature will be used as a "locus" in the phylogeny.
- sequencing reads

outputs:
- a supermatrix alignment fasta
- a RAxML-style "model" file that contains the names and coordinates of individual loci as assembled into the supermatrix
- basic statistics of the supermatrix

The outputs may be used to infer a supermatrix (concatenation) species tree, or separate gene trees as input for a supertree tool.

## -1. setup

```
conda create --name phylo_pipeline snakemake r-geiger bwa samtools bcftools vcftools raxml bedtools seqtk tabix gffread scipy vcflib openjdk iqtree newick_utils -y
conda activate phylo_pipeline 


git clone https://github.com/mscharmann/phylogeny_pipeline

```

## 0. prepare
- adjust parameters and input files in config.yaml
- set up the input read files to samples map in "data/samples_units_readfiles.txt"

## 1. run pipeline
```
snakemake -j 500 --cluster-config cluster.curnagl.json --cluster "sbatch -p {cluster.partition} -t {cluster.time} -c {cluster.CPUs} --mem={cluster.RAM_memory}" --restart-times 3 --keep-going --rerun-incomplete
```

## 2. simple divergence report for four-fold degenerate and all sites:
```
python scripts/count_distance_at_fourfold_degenerate_sites.py results_processed/supermatrix.fasta
```

## 3. gene tree inference
script will calculate branch supports for the gene trees / partitions using IQtree's Ultra-Fast Bootstrap 1000 ("--ufboot 1000")
```
python scripts/iqtree_for_each_partition_multithreading.py results_processed/supermatrix.fasta results_processed/supermatrix.model.txt gene_trees 48
cat gene_trees/*.tre > results_processed/gene_trees.txt
rm -r gene_trees
```

## 4. Recommended: collapse short branches into polytomies 
RAxML will represent polytomies in NEWICK format as sequential dichotomies with very short branch lengths (1e-6). However, the super-tree tool ASTRAL will only use the tree topology, and will be misled by RAxML's representation of polytomies. Thus, to re-code the polytomies as actual NEWICK polytomies, use r-gieger:

```
Rscript scripts/Rscript_collapse_polytomies.txt results_processed/gene_trees.txt results_processed/gene_trees.collapsed_polytomies.txt


```
Newick utilities to collapse nodes with low support, here 70%, into polytomies. This would most likely also collapse the extremely short branches into polytomies, so the above Rscript (which is slow) may not be necessary:

```
nw_ed results_processed/gene_trees.collapsed_polytomies.txt 'i & b<70' o > results_processed/gene_trees.collapsed_polytomies.collapsed_low_support.txt

```


## 5. Supertree inference
Newer version of ASTRAL by default calculates "only the posterior probability for the main resolution". To get Quartet Support instead, apply option "-t 1"

```
wget https://github.com/smirarab/ASTRAL/raw/master/Astral.5.7.8.zip
unzip Astral.5.7.8.zip

java -jar Astral/astral.5.7.8.jar -t 1 -i results_processed/gene_trees.collapsed_polytomies.collapsed_low_support.txt -o results_processed/ASTRAL.tre
```

## 6. Supermatrix (concatenation) species tree
```
raxmlHPC-PTHREADS -T 24 -p 12345 -s results_processed/supermatrix.fasta -n XFU -m GTRCAT -f J
cat RAxML_fastTreeSH_Support.XFU | perl -p -i -e 's/:(\d+\.\d+)\[(\d+)\]/$2:$1/g' >	results_processed/supermatrix_tree.RAxML_fastTreeSH_Support.tre		
rm *XFU*				
```