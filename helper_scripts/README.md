mergeids.py: concatenate two sequences in a fasta based on the sequence header

readlengths.py: print out the length of each sequence in a fasta file

removeids.py: remove all sequences with the IDs given by the user

sortbysize.py: sort the sequences from large to small (assumes all sequences have different lengths so only suitable for files where this is the case)

takelongest.py: take the longest N sequences from a fasta file (assumes the sequences in this range of largest sequences have different lengths)

removelengths.py: remove all reads above or below a certain length

takelongestchrom+extra.py: take the longest N sequences from a fasta file, labels them as chromosomes and the also keeps the remaining n largest sequences (assumes the sequences in this range of largest sequences have different lengths)

uniquecount.py: counts the number of unique values in the first collumn of a tsv file

For information on running the scripts use the -h while running them. All scripts print their output to the terminal.
