from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser(
                    prog='removelengths',
                    description='Script that removes any sequence over a certain length from the fasta',
                    epilog='')

parser.add_argument('input',help="path to the fasta file")
parser.add_argument('threshold',help="size threshold, any sequence longer than this will be filtered out")
parser.add_argument('-i','--invert', action='store_true',help="remove sequences shorten than the threshold as opossed to larger than the threshold")
args = parser.parse_args()
inputfile=args.input
threshold=int(args.threshold)
i=args.invert

if i:
    with open(inputfile) as fasta_file:
        for record in SeqIO.parse(fasta_file, "fasta"):
            if len(record.seq) > threshold:
                print(f">{record.id}\n{record.seq.replace('*','')}")
else:
    with open(inputfile) as fasta_file:
        for record in SeqIO.parse(fasta_file, "fasta"):
            if len(record.seq) < threshold:
                print(f">{record.id}\n{record.seq.replace('*','')}")


