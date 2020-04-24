# TO DO:
# clean up

import sys
import itertools
import argparse
import os
import gzip 

def get_args():
        parser = argparse.ArgumentParser(description='Trim t000_.fasta and 3-mer and 9-mer fragments')
        parser.add_argument("--seq", type=str, default="t000_full.fasta", dest="fseq", help="full sequence [default: %default]")
        parser.add_argument("--3mers", type=str, default="aat000_03_05.200_v1_3_full.txt.gz", dest="f3mers", help="full length 3mers [default: %default]")
        parser.add_argument("--9mers", type=str, default="aat000_09_05.200_v1_3_full.txt.gz", dest="f9mers", help="full length 9mers [default: %default]")
        parser.add_argument("--start", type=int, default="0", dest="istart", help="first amino acid in trimmed sequence [default: %default]")
        parser.add_argument("--stop", type=int, default="0", dest="istop", help="last amino acid in trimmed sequence [default: %default]")
        args = parser.parse_args()

        fseq    = args.fseq
        f3mers  = args.f3mers
        f9mers  = args.f9mers
        istart  = args.istart
        istop   = args.istop

        return fseq, f3mers, f9mers, istart, istop

def get_file_len(fname):
    with gzip.open(fname, "rb") as f:
        for i, l in enumerate(f):
            pass
    return i + 1

### MAIN ###

# Get command line arguments
fseq, f3mers, f9mers, istart, istop = get_args()

# set trimmed filenames
tseq   = "t000_.fasta"
t3mers = "aat000_03_05.200_v1_3.txt.gz"
t9mers = "aat000_09_05.200_v1_3.txt.gz" 

# Remove any existing trimmed files
os.remove(tseq)   if os.path.exists(tseq) else None
os.remove(t3mers) if os.path.exists(t3mers) else None
os.remove(t9mers) if os.path.exists(t9mers) else None

### Trim t000.fasta ###

with open(fseq, "r") as ffull:
    with open(tseq, "a") as ftrimmed:
        for line in ffull:
            if (">") in line:
                line = line.rstrip()
            else:
                line = line.rstrip()
                line = line[istart-1:istop]
            ftrimmed.write(line + '\n')

ffull.close()
ftrimmed.close()

### Trim 3-mers ###

lookup = 'position'

# Generate a list of the start positions for each residue section in the 3-mers.
start = [] 
with gzip.open(f3mers, "rb") as myFile:
    for num, line in enumerate(myFile, 1):
        if lookup in line:
            start.append(num)

# Generate a list of the end positions for each residue section.
stop = []
iterStart = iter(start)
next(iterStart)
for item in iterStart:
    stop.append(int(item) - 1)

flen = get_file_len(f3mers)
stop.append(flen)

icount = 1
with gzip.open(f3mers, "rb") as ffull:
    with gzip.open(t3mers, "ab") as ftrimmed:
        for line in itertools.islice(ffull, start[istart - 1] - 1, stop[istop - 3]):
            if ("position") in line:
                line = ' position:           '+ str(icount) + ' neighbors:          200 '
                icount += 1
            else:
                line = line.rstrip("\n")

            ftrimmed.write(line + '\n')

ffull.close()
ftrimmed.close()

### Trim 9-mers ###

lookup = 'position'

# Generate a list of the start positions for each residue section in the 9-mers.
start = []
with gzip.open(f9mers, "rb") as myFile:
    for num, line in enumerate(myFile, 1):
        if lookup in line:
            start.append(num)

# Generate a list of the end positions for each residue section.
stop = []
iterStart = iter(start)
next(iterStart)
for item in iterStart:
    stop.append(int(item) - 1)

flen = get_file_len(f9mers)
stop.append(flen)

icount = 1
with gzip.open(f9mers, "rb") as ffull:
    with gzip.open(t9mers, "ab") as ftrimmed:
        for line in itertools.islice(ffull, start[istart - 1] - 1, stop[istop - 9]):
            if ("position") in line:
                line = ' position:           '+ str(icount) + ' neighbors:          200 '
                icount += 1
            else:
                line = line.rstrip("\n")

            ftrimmed.write(line + '\n')

ffull.close()
ftrimmed.close()       
