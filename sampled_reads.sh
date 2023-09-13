#!/bin/bash


sample=$1 # data path for bam files
ratio=$2 # sampled fraction
output=$3 # output path

prefix=`basename $sample | sed 's/.bam//g'`

samtools view -@ 8 -s $ratio -b  -o ${output}/${prefix}.${ratio}.bam  $sample
samtools index -@ 8 ${output}/${prefix}.${ratio}.bam
bamCoverage --bam ${output}/${prefix}.${ratio}.bam -o ${output}/${prefix}.${ratio}.bw --binSize 10 --normalizeUsing CPM  --effectiveGenomeSize 2652783500  --numberOfProcessors 8
