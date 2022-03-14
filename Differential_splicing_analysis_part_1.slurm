#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --time=0-2:00:00


module load foss/2019b
module load python/2.7.16
module load samtools/1.9
module load bedtools/2.27.1
module load rmats/4.1.0


echo "AL08_hg19_VaryRL"
echo " "

data="/data/projects/punim1320/data/"
result="/data/projects/punim1320/results/"

gtf=$data$"Homo_sapiens.GRCh37.87.gtf"

CDK4i_72hr=$data$"PD991_72hr.txt"
PRMT5i_72hr=$data$"CTX_72hr.txt"
DMSO_72hr=$data$"DMSO_72hr.txt"
CDK4i_6d=$data$"PD991_6d.txt"


AL=8

cd $result
rm -rf RL149_AL08_hg19_VaryRL
rm -rf RL149_AL08_hg19_VaryRL
mkdir RL149_AL08_hg19_VaryRL
cd $result$"RL149_AL08_hg19_VaryRL"



Compare1="PRMT5i_72hr-DMSO_72hr"
echo $Compare1

rm -rf PRMT5i_72hr_DMSO_72hr

echo "START"
RNASeq-MATS.py --gtf $gtf --b1 $PRMT5i_72hr --b2 $DMSO_72hr --tmp PRMT5i_72hr_DMSO_72hr --od PRMT5i_72hr_DMSO_72hr -t paired --readLength 149 --variable-read-length --anchorLength $AL --nthread 4 --cstat 0.0001 
echo "DONE"
echo " "


Compare2="CDK4i_72hr-DMSO_72hr"
echo $Compare2

rm -rf CDK4i_72hr_DMSO_72hr

echo "START"
RNASeq-MATS.py --gtf $gtf --b1 $CDK4i_72hr --b2 $DMSO_72hr --tmp CDK4i_72hr_DMSO_72hr --od CDK4i_72hr_DMSO_72hr -t paired --readLength 149 --variable-read-length --anchorLength $AL --nthread 4 --cstat 0.0001
echo "DONE"
echo " "


Compare3="CDK4i_6d-DMSO_72hr"
echo $Compare3

rm -rf CDK4i_6d_DMSO_72hr

echo "START"
RNASeq-MATS.py --gtf $gtf --b1 $CDK4i_6d --b2 $DMSO_72hr --tmp CDK4i_6d_DMSO_72hr --od CDK4i_6d_DMSO_72hr -t paired --readLength 149 --variable-read-length --anchorLength $AL --nthread 4 --cstat 0.0001
echo "DONE"
echo " "


Compare4="CDK4i_72hr-PRMT5i_72hr"
echo $Compare4

rm -rf CDK4i_72hr_PRMT5i_72hr

echo "START"
RNASeq-MATS.py --gtf $gtf --b1 $CDK4i_72hr --b2 $PRMT5i_72hr --tmp CDK4i_72hr_PRMT5i_72hr --od CDK4i_72hr_PRMT5i_72hr -t paired --readLength 149 --variable-read-length --anchorLength $AL --nthread 4 --cstat 0.0001
echo "DONE"
echo " "


Compare5="CDK4i_6d-PRMT5i_72hr"
echo $Compare5

rm -rf CDK4i_6d_PRMT5i_72hr

echo "START"
RNASeq-MATS.py --gtf $gtf --b1 $CDK4i_6d --b2 $PRMT5i_72hr --tmp CDK4i_6d_PRMT5i_72hr --od CDK4i_6d_PRMT5i_72hr -t paired --readLength 149 --variable-read-length --anchorLength $AL --nthread 4 --cstat 0.0001
echo "DONE"
echo " "


Compare6="CDK4i_6d-CDK4i_72hr"
echo $Compare6

rm -rf CDK4i_6d_CDK4i_72hr

echo "START"
RNASeq-MATS.py --gtf $gtf --b1 $CDK4i_6d --b2 $CDK4i_72hr --tmp CDK4i_6d_CDK4i_72hr --od CDK4i_6d_CDK4i_72hr -t paired --readLength 149 --variable-read-length --anchorLength $AL --nthread 4 --cstat 0.0001
echo "DONE"
echo " "




echo " "
echo "All completed"



