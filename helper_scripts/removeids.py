from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser(
                    prog='removeids',
                    description='Script that gives removes a set of sequences from a fasta file based on the IDs',
                    epilog='')

parser.add_argument('input',help="path to the fasta file")
parser.add_argument('ids',help="list of ID's to ignore (format: ID1 ID2 ID3 etc.)", nargs="*")

args = parser.parse_args()

inputfile=args.input
idlist=args.ids

with open(inputfile) as fasta_file:
    for record in SeqIO.parse(fasta_file, "fasta"):
        if record.id not in idlist:
            print(f">{record.id}\n{record.seq}")
