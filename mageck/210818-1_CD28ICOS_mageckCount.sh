#!/bin/bash                         #-- what is the language of this shell
#                                  #-- Any line that starts with #$ is an instruction to SGE
#$ -S /bin/bash                     #-- the shell for the job
#$ -o /wynton/group/ye/cody/logs/                        #-- output directory (fill in)
#$ -e /wynton/group/ye/cody/logs/                        #-- error directory (fill in)
#$ -cwd                            #-- tell the job that it should start in your working directory
#$ -r y                            #-- tell the system that if a job crashes, it should be restarted
#$ -j y                            #-- tell the system that the STDERR and STDOUT should be joined
#$ -l mem_free=15G                  #-- submits on nodes with enough free memory (required)
#$ -l scratch=5G         #-- SGE resources (home and scratch disks)
#$ -l h_rt=12:00:00                #-- runtime limit (see above; this requests 24 hours)
#$ -m ea                           #--email when done
#$ -M cody.mowery@gladstone.ucsf.edu        #--email
#$ -pe smp 1 
#$ -t 1-4 

ulimit -n 16000

conda activate screens

script_path=/wynton/group/ye/cody/screens/git_scripts/tadScreen/finalScripts/mageck/
donor=$(awk 'NR==i{print $1}' i=${SGE_TASK_ID} ${script_path}210818-1_CD28ICOS_mageckSamplesheet.txt)
target=$(awk 'NR==i{print $2}' i=${SGE_TASK_ID} ${script_path}210818-1_CD28ICOS_mageckSamplesheet.txt)
lo=$(awk 'NR==i{print $3}' i=${SGE_TASK_ID} ${script_path}210818-1_CD28ICOS_mageckSamplesheet.txt)
hi=$(awk 'NR==i{print $4}' i=${SGE_TASK_ID} ${script_path}210818-1_CD28ICOS_mageckSamplesheet.txt)
echo ${donor} ${target}

ref=${script_path}UniqueLibraryRef.csv

fq_path=/wynton/group/ye/cody/screens/data/210818-1/210818-1_CD28_ICOS/fastqs/
outdir=/wynton/group/ye/cody/screens/data/210818-1/210818-1_CD28_ICOS/mageckFinal/

cd ${outdir}

mageck count -l ${ref} \
 -n ${donor}_${target} \
 --sample-label ${donor}_${target}_LO,${donor}_${target}_HI \
 --fastq ${fq_path}*${lo}_*R1_001.fastq.gz ${fq_path}*${hi}_*R1_001.fastq.gz

cp ${outdir}${donor}_${target}*count.txt ${script_path}/210818-1_mageckOuts/

qstat -j $JOB_ID            #-- helpful run info for debugging etc
