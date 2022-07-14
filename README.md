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
conda create --name phylo_pipeline snakemake r-geiger bwa samtools bcftools vcftools raxml bedtools seqtk tabix gffread scipy vcflib -y
conda activate phylo_pipeline 


git clone https://github.com/mscharmann/phylogeny_pipeline

```

## 0. prepare
- adjust parameters and input files in config.yaml
- set up the input read files to samples map in "data/samples_units_readfiles.txt"

## 1. run pipeline
```
snakemake -j 500 --cluster-config cluster.axiom.json --cluster "sbatch -p {cluster.partition} -t {cluster.time} -c {cluster.CPUs} --mem={cluster.RAM_memory}" --restart-times 3 --keep-going --rerun-incomplete
```

## 2. simple divergence report for four-fold degenerate and all sites:
```
python scripts/count_distance_at_fourfold_degenerate_sites.py results_processed/supermatrix.fasta
```

## 3. gene tree inference
script will NOT calculate branch supports for the gene trees / partitions.
```
python scripts/raxml_for_each_partition_multithreading.py results_processed/supermatrix.fasta results_processed/supermatrix.model.txt gene_trees 48
cat gene_trees/*.tre > gene_trees.txt
rm -r gene_trees
```

## 4. Recommended: collapse short branches into polytomies 
RAxML will represent polytomies in NEWICK format as sequential dichotomies with very short branch lengths (1e-6). However, the super-tree tool ASTRAL will only use the tree topology, and will be misled by RAxML's representation of polytomies. Thus, to re-code the polytomies as actual NEWICK polytomies, use r-gieger:

```
Rscript scripts/Rscript_collapse_polytomies.txt gene_trees.txt > gene_trees_for_ASTRAL.txt
```

## 5. Supertree inference

```
wget https://github.com/smirarab/ASTRAL/raw/master/Astral.5.7.8.zip
unzip Astral.5.7.8.zip

java -jar Astral/astral.5.7.8.jar -i gene_trees_for_ASTRAL.txt -o ASTRAL.tre
```

## 6. Supermatrix (concatenation) species tree
```
raxmlHPC-PTHREADS -T 24 -p 12345 -s results_processed/supermatrix.fasta -n XFU -m GTRCAT -f J
cat RAxML_fastTreeSH_Support.XFU | perl -p -i -e 's/:(\d+\.\d+)\[(\d+)\]/$2:$1/g' >	supermatrix_tree.RAxML_fastTreeSH_Support.tre		
rm *XFU*				
```