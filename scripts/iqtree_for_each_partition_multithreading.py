# iqtree_for_each_partition_multithreading.py
# Mathias Scharmann
# 2022-07



"""
Input:
- aligment in FASTA format
- raxml partition format file

Output: 
- a directory with trees estimated separately for each partition!


module load gdc
module load python/2.7
PATH=$PATH:/cluster/project/gdc/people/schamath/tools


"""

import sys, os
import subprocess

###

def is_non_zero_file(fpath):  
    return os.path.isfile(fpath) and os.path.getsize(fpath) > 0

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
					seq += line.strip("\n")
	
		# append last seq
		outlines.append(seq)
	
	i=0
	j=1
	out_dict = {}
	for x in range(int(len(outlines)/2.0)):
		out_dict[outlines[i]] = outlines[j]
		i += 2
		j += 2
	
	return out_dict


def write_phylip_clean (cl_dict, outfile):
		
	with open(outfile, "w") as OUTFILE:
		ntaxa = len(cl_dict.keys())
		len_align = len(list(cl_dict.values())[0]) # gets value of first element in dictionary -> this OK since all seqs have same length
		header = str(ntaxa) + " " + str(len_align) + "\n"
		OUTFILE.write(header)
		for sample, seq in cl_dict.items():
			out_line = sample + "    " + seq + "\n"
			OUTFILE.write(out_line)
	OUTFILE.close()


def split_alignment_and_do_iqtree ( alignment_file, partitions_file, outdir, parallel_processes ):
	
	master_aln = read_wrapped_or_unwrapped_fasta( alignment_file )
	
	partitions_dict = {}
	with open(partitions_file, "r") as INFILE:
		for line in INFILE:
			name = line.split(",")[1].split("=")[0]	
			start_idx = int(line.split("=")[1].split("-")[0]) -1 # zero-based offset for python idxes
			end_idx = int(line.split("=")[1].split("-")[1] ) #+ 1 # string/list slicing in python: right end idx NOT included in slice
			partitions_dict[name] = [start_idx, end_idx]								
	
	cnt = 0
	processes = set()
	max_processes = parallel_processes
	num_cores_per_process = num_cores = 1
	for name, idxes in partitions_dict.items():
		cnt += 1		
		outdict = {tax: master_aln[tax][partitions_dict[name][0] : partitions_dict[name][1]] for tax in master_aln.keys() }		
		# drop taxa that have only undetermined characters "-", otherwise RAxML will complain
		clean_dict = {k:v for k,v in outdict.items() if set(v) != set("-") }
		
		# sometimes, alignment is slightly incorrect and contains ONLY missing data for a locus; catch these errors and continue. Also, gene trees aren't informative if there is not at least 4 tips on them.
		if len(clean_dict.keys()) >= 4:
			
			# also, don't run again if there is already a treefile with content in it:
			if not is_non_zero_file("./" + outdir + "/" + name + ".raxml.tre"):

				write_phylip_clean (clean_dict , outdir + "/" +  name + ".phy")
				
				run_name = "run_" + str(cnt) + "_"
				
				cmd = " ".join ( ["cd", outdir,";", "iqtree","-T",str(num_cores),"-s",\
								name + ".phy" ,"--quiet", "--ufboot 1000" , ";",\
								"mv", name + ".phy.treefile", name + ".iqtree.tre", ";",\
								"rm", name + ".phy*"] )
	
	#			print cmd
				processes.add(subprocess.Popen([cmd, run_name], shell = True))
				print(run_name)
				if len(processes) >= max_processes:
					os.wait()
					processes.difference_update([ p for p in processes if p.poll() is not None])

########## main

if len(sys.argv) != 5:
	print ("""usage: python iqtree_for_each_partition_multithreading.py alignment_file partitions_file outdir total_cores\n
IQTREE will be called on 1 core for each parallel job""")
	
	exit()

alignment_file = sys.argv[1]
partitions_file = sys.argv[2]
outdir = sys.argv[3]
try:
	os.mkdir( outdir )
except OSError:
	None

total_cores = int( sys.argv[4] )
num_cores_per_process = 1
parallel_processes = total_cores / num_cores_per_process

split_alignment_and_do_iqtree ( alignment_file, partitions_file, outdir, parallel_processes )


