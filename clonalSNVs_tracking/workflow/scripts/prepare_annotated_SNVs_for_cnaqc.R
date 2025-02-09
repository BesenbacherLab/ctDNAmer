sink(snakemake@log[[1]], append=TRUE) # logging

# packages
library(tidyverse)
library(stringr)

# helper functions
find_header = function(file) {
    header_pref_vec = c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT")
    data = file(file, "r")
    while ( TRUE ) {
        line = readLines(data, n = 1)
        if ( length(line) == 0 ) {
            break
        }
        if (identical(str_split_fixed(line, "[ \t]", 11)[1:9], header_pref_vec)){
            h = str_split_fixed(line, "[ \t]", 11)[1:11]
            if (grepl("buffycoat", h[10], fixed = TRUE) == TRUE){
                h[10:11] = c("buffycoat", "tumor")
            } else{
                h[10:11] = c("tumor", "buffycoat")
            }
            h[1] = "CHROM"
            return(h)
        }
      }
    close(data)
}

get_metadata <- function(data){
    impact <- str_split_fixed(str_split_fixed(str_split_fixed(data$INFO, ";", 15)[, 15], "=", 2)[, 2], "[|]", 10)[, 6]
    gene_label <- str_split_fixed(str_split_fixed(str_split_fixed(data$INFO, ";", 15)[, 15], "=", 2)[, 2], "[|]", 10)[, 9]
    n_alt_alleles <- ifelse(substr(data$FORMAT, 1, 11) == "GT:AD:AF:DP", str_split_fixed(str_split_fixed(data$tumor, ":", 5)[, 2], ",", 2)[, 2], NA)
    depth <- ifelse(substr(data$FORMAT, 1, 11) == "GT:AD:AF:DP", str_split_fixed(data$tumor, ":", 5)[, 4], NA)
    return(tibble(impact = impact, 
                  driver_label = gene_label, 
                  NV = as.numeric(n_alt_alleles), 
                  DP = as.numeric(depth)))
}

# input files
mutect_ann_file = snakemake@input[["m_calls_annot"]]
print(paste0("Input filename (annotated mutect calls): ", mutect_ann_file))

mutect_ann_coding_file = snakemake@input[["m_calls_annot_coding"]]
print(paste0("Input filename (annotated mutect calls - coding only): ", mutect_ann_coding_file))

# output files
output_file_no_ann <- snakemake@output[["m_calls_cnaqc"]]
print(paste0("Output file (no annotation): ", output_file_no_ann))
output_file_ann <- snakemake@output[["m_calls_annot_cnaqc"]]
print(paste0("Output file (annotated): ", output_file_ann))

# get Mutect2 calls header
header = find_header(mutect_ann_file)
print("SNV calls header")
print(header)

# process all SNV annotations
all_ann <- read.table(mutect_ann_file, stringsAsFactors=FALSE, comment.char="#")
colnames(all_ann) <- header

all_ann <- cbind(all_ann, get_metadata(all_ann))
all_ann <- all_ann %>% rowwise() %>% mutate(VAF = NV/DP, from = POS - 1) %>% select(CHROM, from, POS, REF, ALT, DP, NV, VAF) %>% 
    rename(to = POS, chr = CHROM, ref = REF, alt = ALT)
print("Head and dim of all mutations annotations")
print(head(all_ann, n = 3))

# dummy checks
print(paste0("Number of mutations: ", nrow(all_ann)))
print(paste0("Number of NA values in total: ", sum(is.na(all_ann))))
print(paste0("Number of NA values in VAF column: ", sum(is.na(all_ann$VAF))))
# write results
write.table(all_ann, output_file_no_ann, row.names=FALSE, col.names=TRUE, sep="\t", quote = FALSE)

# process coding SNV annotations
coding_ann <- read.table(mutect_ann_coding_file, stringsAsFactors=FALSE,comment.char="#")
colnames(coding_ann) <- header
coding_ann <- coding_ann %>% filter(grepl('CSQ', INFO))
coding_ann <- cbind(coding_ann, get_metadata(coding_ann))
coding_ann <- coding_ann %>% rowwise() %>% mutate(VAF = NV/DP, from = POS - 1) %>% select(CHROM, from, POS, REF, ALT, DP, NV, VAF, impact, driver_label) %>%
    rename(chr = CHROM, to = POS, ref = REF, alt = ALT)
coding_ann <- coding_ann %>% filter(impact %in% c("HIGH"))
print("Head and dim of coding mutations annotations")
print(head(coding_ann, n = 3))
print(dim(coding_ann))

all_ann <- left_join(all_ann, coding_ann, by = c("chr", "from", "to", "ref", "alt", "DP", "NV", "VAF"))
all_ann <- all_ann %>% mutate(is_driver = ifelse(is.na(impact), FALSE, TRUE)) %>% select(-impact) %>% select(chr, from, to, ref, alt, DP, NV, VAF, is_driver, driver_label)

# dummy checks
print(paste0("Number of mutations: ", nrow(all_ann)))
print(paste0("Number of non-NA values in driver label: ", sum(!(is.na(all_ann$driver_label)))))

print("Head and dim of all mutations with high-impact coding mutations annotated")
print(head(all_ann, n = 3))
# write results
write.table(all_ann, output_file_ann, row.names=FALSE, col.names=TRUE, sep="\t", quote = FALSE)

sink()