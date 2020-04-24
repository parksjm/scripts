#!/bin/csh -f

if ($#argv == 1) then
  set job = $argv[1]
  echo "getting robetta output for job # $job"
else
  echo "Usage: [script] [job #]"
endif

curl "http://www.robetta.org/downloads/fragments/$job/t000_.fasta" > t000_full.fasta
curl "http://www.robetta.org/downloads/fragments/$job/aat000_03_05.200_v1_3" > aat000_03_05.200_v1_3_full.txt
curl "http://www.robetta.org/downloads/fragments/$job/aat000_09_05.200_v1_3" > aat000_09_05.200_v1_3_full.txt

gzip aat000_03_05.200_v1_3_full.txt
gzip aat000_09_05.200_v1_3_full.txt

# Change multi line fasta to single line fasta
awk '/^>/{print s? s"\n"$0:$0;s="";next}{s=s sprintf("%s",$0)}END{if(s)print s}' t000_full.fasta > tmp
mv tmp t000_full.fasta 
