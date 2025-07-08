## Clear the environmental space here to start affresh 
## Clear environmental space
rm(list = ls())

## Set the directory with the expression datasets here 
## Set path to the staroutfiles here
setwd('/Volumes/Disk/PhD/Data/PhD_data/Brain/GZ11_star_output')

# This is the output folder for the final analysis
output_folder = '/Volumes/Disk/PhD/Data/PhD_data/March_03_25_Final_Analysis/DGE_Files/'
output_dir = '/Users/gwk/Desktop/PhD/Data/PhD_data/March_03_25_Final_Analysis/normalised_data_sets'

## Source user defined files containing some useful functions here
source("/Users/gwk/Desktop/Bioinformatics/bulk-rnaseq-zebrafish-mouse-human/scripts/functionsneededforanalysis.R")

## Read the mapped and quantified files into r for further analysis here 
file1 <- "1_S1_ReadsPerGene.out.tab"
dr <- '/Volumes/Disk/PhD/Data/PhD_data/Brain/GZ11_star_output/star'

## load the data into R for further analysis here 
# This is output from star - these are files with ReadsPerGene.out.tab
countMatrix <- staroutput_preprocessing(dr, file1)

## Read in the data from the sample informatio here 
# If no sample information exists, you can create to match the countMatrix
samples <- read.csv("sample_information.csv", row.names = 1)

## Explore the loaded data here
head(countMatrix)
head(samples, 15)

## Fix the column names in the count matrix to match the sample information rownames
row_col_names <- as.vector(paste0("sample",rownames(samples)))
rownames(samples) <- row_col_names
colnames(countMatrix) <- row_col_names

## Makesure that sample rownames match the countmatrix colnames here
# The lines below must return true for all
# if false, please investigate the rownames match the countMatrix colnames
all(row.names(samples) %in% colnames(countMatrix))
all(row.names(samples) == colnames(countMatrix))

# Write data to a seperate file so that you do not have to run this step anymore
