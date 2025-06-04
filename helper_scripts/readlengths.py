from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser(
                    prog='firstread',
                    description='Script that gives read lengths from a fasta file',
                    epilog='')

parser.add_argument('input',help="path to the fasta file")
args = parser.parse_args()
inputfile=args.input

with open(inputfile) as fasta_file:
    for record in SeqIO.parse(fasta_file, "fasta"):
        print(f"{record.id}\t{len(record.seq)}")

