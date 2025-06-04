from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser(
                    prog='mergeids',
                    description='Script that merges 2 given ids from a fasta file',
                    epilog='')

parser.add_argument('input',help="path to the fasta file")
parser.add_argument('id1',help="the first header id (don't include the >)")
parser.add_argument('id2',help="the second header (don't include the >)")
args = parser.parse_args()
inputfile=args.input
id1=args.id1
id2=args.id2
mixhead=f"{id1}+{id2}"
mixseq=""
outtext=""
with open(inputfile) as fasta_file:
    for record in SeqIO.parse(fasta_file, "fasta"):
        if record.id==id1:
            mixseq+=record.seq
        elif record.id==id2:
            mixseq+=record.seq
        else:
            outtext+=f">{record.id}\n{record.seq}\n"

outtext+=f">{mixhead}\n{mixseq}\n"

print(outtext)



