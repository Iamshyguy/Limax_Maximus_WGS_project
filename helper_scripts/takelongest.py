from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser(
                    prog='longest',
                    description='Script that takes the longest n sequences from a fasta file',
                    epilog='')

parser.add_argument('input',help="path to the fasta file")
parser.add_argument('n',help="amount of sequences to take")

args = parser.parse_args()
inputfile=args.input
n=int(args.n)
lengths=[]

outtext=""
with open(inputfile) as fasta_file:
    for record in SeqIO.parse(fasta_file, "fasta"):
        lengths.append(len(record.seq))
    lengths.sort()
    num=n*-1
    keep=lengths[num:]

with open(inputfile) as fasta_file:
    for record in SeqIO.parse(fasta_file, "fasta"):
        if len(record.seq) in keep:
            outtext+=f">{record.id}\n{record.seq}\n"
print(outtext)
