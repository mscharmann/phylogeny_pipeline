import sys

thresh = float(sys.argv[1])
INFASTA = sys.argv[2]

#thresh = 0.7
#INFASTA = "window_1.fasta"

inlist = []
with open(INFASTA, "r") as F:
	for line in F:
		if line.startswith(">"):
			seqid = line.strip("\n").strip(">")
		else:
			seq = line.strip("\n")
			inlist.append([seqid, seq])

		
ntax = len(inlist)
min_tax = int(thresh*ntax) 

seqids = [s[0] for s in inlist]

outlist = [ [id,""] for id in seqids ]
for col in range(len(inlist[0][1])):
	chars = [x[col] for x in [s[1] for s in inlist] if not x[col] == "-"]
	if len(chars) >= min_tax:
		for idx in range(len(seqids)):
			outlist[idx][1] += inlist[idx][1][col]

for s in outlist:
	sys.stdout.write(">" + s[0] + "\n" + s[1] + "\n")
