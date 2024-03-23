#!/bin/bash -l
#$ -N pipe4C_pipeline
#$ -cwd
#$ -pe smp 8 -R y -binding linear:8
#$ -l h_vmem=12G
#$ -e pipe4C_pipeline.log
#$ -o pipe4C_pipeline.log
#$ -j yes
#$ -l h_rt=96:00:00


conda activate pipe4C
use Samtools
use Bowtie2
use .r-3.6.0-bioconductor

Rscript /seq/epiprod02/sgaldon/Zeyu/4C/pipe4C/pipe4C.R --vpFile=./config_files/VPinfo.txt --fqFolder=./fastqs/ --outFolder=./pipe4Cout --cores 8 --plot --wig
