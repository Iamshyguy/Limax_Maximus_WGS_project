from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser(
                    prog='countunique',
                    description='Script that counts the unique values in the nth column of a tsv (0 indexed!)',
                    epilog='')

parser.add_argument('input',help="path to the tsv file")
parser.add_argument('n',help="the column to count (0 indexed!)")
args = parser.parse_args()
inputfile=args.input
n=int(args.n)
count=0
listo=[]
with open(inputfile,"r") as file:
    lines=file.readlines()
    for line in lines:
        id=line.split("\t")[n]
        if id not in listo:
            count+=1
            listo.append(id)
print(count)




