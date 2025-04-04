#!/usr/bin/bash -l

## State the number of cores 
#$ -pe mpi 5

## Request for the number of threads
#$ -pe smp 8

## State the amount of RAM per core
#$ -l mem=10G

## Select the run time
#$ -l h_rt=48:00:00

## Request 50 gigabyte of TMPDIR 
#$ -l tmpfs=30G

## Specify the working directory here
#$ -wd /home/uczrgwk/Scratch/mouse/reference

## Load the modules that are required for teh alignment here 
## Unload the gcc-libs module before attempting the following 
module unload gcc-libs

## Load the modules that are required for use to download sra reads
module load gcc-libs/10.2.0
module load sra-tools/3.0.6/gnu-10.2.0
module load star/2.7.10b/sbindist

##Activate the conda environment here pointing to the tools needed
conda activate /home/uczrgwk/miniconda3/envs/BulkRna

## The script to be submitted goes here