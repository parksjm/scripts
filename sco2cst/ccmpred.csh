#!/bin/csh -f 

#set msa = 'PF16998'
set msa = '278'
echo 'msa = ' $msa

python convert_alignment.py $msa.fas fasta $msa.aln

#/home/jpc/bin/ccmpred -t 32 -n 100 -e 0 -A $msa.aln $msa.mtx > ccmpred.out
#/Users/jpc/bin/ccmpred -n 100 -e 0 -A $msa.aln $msa.mtx > ccmpred.out
echo 'skipping ccmpred'


python top_couplings.py $msa.mtx -s 3 -n 500 > scores.out
