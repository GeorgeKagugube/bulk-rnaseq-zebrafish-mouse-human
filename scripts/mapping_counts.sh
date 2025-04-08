## Set up a shebang here 
#!/shared/ucl/apps/python/bundles/python39-6.0.0/venv/bin/python

## Import the required modules
import subprocess as sp
import os 
import sys

## Set up the directories for the accession numbers
ctl_dir = '/home/uczrgwk/Scratch/mouse/cp/ctl'
mut_dir = '/home/uczrgwk/Scratch/mouse/cp/mut'
outputdir = '/home/uczrgwk/Scratch/mouse/cp/results'
ref = "/home/uczrgwk/Scratch/mouse/reference"

## Download another reference genom ehere 
wget -P $ref https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.27_GRCm39/GCF_000001635.27_GRCm39_genomic.fna.gz
wget -P $ref https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.27_GRCm39/GCF_000001635.27_GRCm39_genomic.gtf.gz

"""
## Run fasterq dump command to download the data required here
fasterq-dump <accession number> \
			 --outdir $outputdir \
			 --split-files

## Run fastqc on the reads here
fastqc <fastq files as a list> -o [output dir] \
	   --threads ${NSLOTS/2}

## Once all fastqc is run, amalgamate the data into one readible file using multiqc
multiqc .
"""

## Create a reference for the orginism of interest
## Run star for the creation of the index here
build_index = """STAR --runMode genomeGenerate \
	 --genomeDir ${ref} \
	 --genomeFastaFiles ${ref}/*.fna.gz \
	 --sjdbGTFfile ${ref}/*.gtf.gz \
	 --runThreadN ${NSLOTS}
	 """

## Run Python call 
sb.run([build_index], shell=True)

## Run the first pass STAR
## run star alignment on the assembly here 
round1 = """STAR --runThreadN ${NSLOTS} \
	 --genomeDir ${ref} \
	 --readFilesIn ${forwardreads} ${reverseread} \
	 --outFileNamePrefix ${outdir1}/${outputfilename}_ \
     --readFilesCommand zcat \
	 --quantMode GeneCounts \
	 --outSAMtype BAM SortedByCoordinate 
	 """

## run star alignment on the assembly here 
## run the second pass now 
round2 = """STAR --runThreadN ${NSLOTS} \
	 --genomeDir $ref \
	 --readFilesIn $forwardreads $reverseread  \
	 --readFilesCommand zcat \
	 --outFileNamePrefix ${outdir2}/${outputfilename}_ \
	 --quantMode GeneCounts \
	 --outSAMtype BAM SortedByCoordinate \
	 --sjdbFileChrStartEnd `find ${outdir1} | grep SJ.out.tab* | sort | tr '\n' ' '`
	 """
