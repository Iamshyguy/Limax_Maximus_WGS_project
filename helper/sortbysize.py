from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser(
                    prog='longest',
                    description='Script sorts sequences by length (assumes lengths are unique, if not, the script wont work)')

parser.add_argument('input',help="path to the fasta file")

args = parser.parse_args()
inputfile=args.input

lengths=[]
lendir={}

outtext=""
with open(inputfile) as fasta_file:
    for record in SeqIO.parse(fasta_file, "fasta"):
        lengths.append(len(record.seq))
        lendir[len(record.seq)]=[record.id,record.seq]

    lengths.sort(reverse=True)
    for length in lengths:
        outtext+=f">{lendir[length][0]}\n{lendir[length][1]}\n"
print(outtext)
