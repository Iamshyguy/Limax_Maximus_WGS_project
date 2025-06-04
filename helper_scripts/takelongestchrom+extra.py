from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser(
                    prog='longest',
                    description='Script that takes the longest n sequences from a fasta file, renaming them to chromosomes. and a number of additional scaffolds',
                    epilog='')

parser.add_argument('input',help="path to the fasta file")
parser.add_argument('n',help="amount of sequences to take")
parser.add_argument('-c','--cutoff',help="amount of additional sequences to take")

args = parser.parse_args()
inputfile=args.input
n=int(args.n)
try:
    cutoff=int(args.cutoff)
except:
    cutoff=0
    print(huh)
else:
    cutoff=int(args.cutoff)
lengths=[]
val=1
outextra=""
outtext=""
with open(inputfile) as fasta_file:
    for record in SeqIO.parse(fasta_file, "fasta"):
        lengths.append(len(record.seq))
    lengths.sort()
    num=n*-1
    if cutoff == 0:
        num2=0
    else:
        num2=cutoff*-1
        num2+=num
    keep=lengths[num:]
    merge=lengths[num2:num]


with open(inputfile) as fasta_file:
    for record in SeqIO.parse(fasta_file, "fasta"):
        if len(record.seq) in keep:
            outtext+=f">chromosome_{val}\n{record.seq}\n"
            val+=1
        if len(record.seq) in merge:
            outextra+=f">{record.id}\n{record.seq}\n"

outtext+=outextra

print(outtext)
