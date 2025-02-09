sink(snakemake@log[[1]], append=TRUE) # logging

# packages
library(sequenza)

# input params
pt_id <- snakemake@params[["pt_id"]] # patient ID
print(paste0("pt_id: ", pt_id))

output_dir <- snakemake@params[["output_dir"]] # output directory
print(paste0("Output directory: ", output_dir))

female <- as.character(snakemake@params[["female"]]) # patinet gender
female_logic <- ifelse(female == "True", TRUE, FALSE)
print(paste0("pt female, logic_var: ", female_logic))

# purity and ploidy limits
cellularity_min <- as.numeric(snakemake@params[["cellularity_min"]])
cellularity_max <- as.numeric(snakemake@params[["cellularity_max"]])
print(paste0("Cellularity limits: ", cellularity_min, ", ", cellularity_max))
ploidy_min <- as.numeric(snakemake@params[["ploidy_min"]])
ploidy_max <- as.numeric(snakemake@params[["ploidy_max"]])
print(paste0("Ploidy limits: ", ploidy_min, ", ", ploidy_max))

# choose chromosomes based on patient gender
if (female_logic){
    print("Patient is female")
    chr_list = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", 
                 "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", 
                 "chr20", "chr21", "chr22", "chrX")
} else {
    chr_list = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", 
                 "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", 
                 "chr20", "chr21", "chr22")
}

input_file = snakemake@input[["seq_binned"]] # read input file
print(paste0("Input filename: ", input_file))

# sequenza.extract: process seqz data, normalization and segmentation
test <- sequenza.extract(input_file, assembly = "hg38", verbose = TRUE, 
                         chromosome.list = chr_list, 
                         gamma = 280, kmin = 300, max.mut.types = 1, min.reads.baf = 5, min.reads = 10)


# sequenza.fit: run grid-search approach to estimate cellularity and ploidy
set.seed(0)
CP <- sequenza.fit(test, female = female_logic, 
                  cellularity = seq(cellularity_min, cellularity_max, 0.01), 
                  ploidy = seq(ploidy_min, ploidy_max, 0.1))

# sequenza.results: write files and plots using suggested or selected solution
sequenza.results(sequenza.extract = test,
    cp.table = CP, sample.id = pt_id,
    out.dir=output_dir, female = female_logic)

sink()