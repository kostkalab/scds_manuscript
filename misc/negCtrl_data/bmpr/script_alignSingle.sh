#!/bin/bash
############################################################
## Generate indexed, sorted alignments using Hisat2 for single
## sample. 

## The code below needs Hisat2 and samtools installed, and Hisat2
## indexes for grcm38 set up. Please modify relevant paths below
## accordingly. 

## Options modified per original article
## https://www.sciencedirect.com/science/article/pii/S2211124717309610?via%3Dihub#app3
############################################################
set -x

############################################################
## arguments
############################################################
srrId=$1
## srrId='SRR3092152'

############################################################
## prelims/project-related
############################################################
PROJ_DIR='./' 
dataSetName='bmpr'

RAWDATA_DIR=$PROJ_DIR/$dataSetName/raw/
PROCDATA_DIR=$PROJ_DIR/$dataSetName/proc/

############################################################
ANNO_DIR='./refData/hisat2/grcm38/' ## modify path before running
refGenHisat2Index=$ANNO_DIR/genome

############################################################
OUTPUT_DIR=$PROCDATA_DIR/$(basename $srrId)/hisat2/
mkdir -p $OUTPUT_DIR

LOGS_DIR=$PROJ_DIR/$dataSetName/logs/ 
mkdir -p $LOGS_DIR

inFile1=$RAWDATA_DIR/${srrId}_pass_1.fastq.gz
inFile2=$RAWDATA_DIR/${srrId}_pass_2.fastq.gz

outputBAM=$OUTPUT_DIR/${srrId}.bam
sortedBAM=$OUTPUT_DIR/${srrId}.sorted.bam

############################################################
## Fastq -> BAM
############################################################
logFile=$LOGS_DIR/"00_align-"${srrId}.out
summaryFile=$OUTPUT_DIR/summary.txt
metricsFile=$OUTPUT_DIR/metrics.txt

date "+Starting file: $srrId %Y-%m-%d %H:%M:%S %n" >$logFile
echo -e "Log file: $logFile " >>$logFile
echo -e "Sample: "$srrId      >>$logFile 

date "+Running HISAT2: %Y-%m-%d %H:%M:%S" >>$logFile
hisat2 \
    -x $refGenHisat2Index \
    -p 4 \
    --no-mixed \
    --no-discordant \
    --no-softclip \
    --summary-file $summaryFile \
    --met-file $metricsFile \
    -1 $inFile1 \
    -2 $inFile2 \
    | samtools view -bS -F 256 -q 20 - > $outputBAM \
    && date "+Finished running HISAT2 %Y-%m-%d %H:%M:%S %n" >>$logFile 2>&1 || exit 1

############################################################
## BAM -> sorted BAM
############################################################
date "+Sort BAM: %Y-%m-%d %H:%M:%S" >>$logFile
samtools sort -@ 2 -O BAM -o $sortedBAM $outputBAM 1>> $logFile 2>&1 || exit 1

############################################################
## sorted BAM -> sorted BAM.bai
############################################################
date "+Index sorted BAM: %Y-%m-%d %H:%M:%S" >>$logFile
samtools index $sortedBAM 1>> $logFile 2>&1 || exit 1

date '+Done Index BAM \%y-\%m-\%d \%H:\%M:\%S' >>$logFile
date "+Done sample: $srrId %Y-%m-%d %H:%M:%S %n" >$OUTPUT_DIR/hisat.done

############################################################
############################################################
############################################################

    








