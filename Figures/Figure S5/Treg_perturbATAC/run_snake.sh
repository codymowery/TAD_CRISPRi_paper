#!/bin/bash                         #-- what is the language of this shell
#                                  #-- Any line that starts with #$ is an instruction to SGE
#$ -S /bin/bash                     #-- the shell for the job
#$ -o /wynton/group/ye/cody/logs/                        #-- output directory (fill in)
#$ -e /wynton/group/ye/cody/logs/                        #-- error directory (fill in)
#$ -cwd                            #-- tell the job that it should start in your working directory
#$ -r y                            #-- tell the system that if a job crashes, it should be restarted
#$ -j y                            #-- tell the system that the STDERR and STDOUT should be joined
#$ -l mem_free=10G                  #-- submits on nodes with enough free memory (required)
#$ -l scratch=10G         #-- SGE resources (home and scratch disks)
#$ -l h_rt=72:00:00                #-- runtime limit (see above; this requests 24 hours)
#$ -m ea                           #--email when done
#$ -M cody.mowery@gladstone.ucsf.edu        #--email
#$ -pe smp 1

ulimit -n 16000

conda activate snakemake_freimerATAC
cd ATAC_Seq_pipeline/

snakemake --default-resources mem_mb=125000 disk_mb=125000 --cores ${NSLOTS} --use-conda --jobs 64 --cluster "qsub -V -b y" --rerun-incomplete


qstat -j $JOB_ID
