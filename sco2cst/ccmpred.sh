#!/bin/bash 

msa='278'
echo 'msa = ' $msa

python convert_alignment.py $msa.fas fasta $msa.aln

#ccmpred=/Users/jpc/bin/ccmpred
#$ccmpred -t 32 -n 100 -e 0 -A $msa.aln $msa.mtx > ccmpred.out

#ccmpred=/home/jpc/bin/ccmpred
#$ccmpred -n 100 -e 0 -A $msa.aln $msa.mtx > ccmpred.out
echo 'skipping ccmpred'


python top_couplings.py $msa.mtx -s 3 -n 500 > scores.out
