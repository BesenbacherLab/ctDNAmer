sink(snakemake@log[[1]], append=TRUE)  # logging

# packages
library(tidyverse)
library(stringr)
library(CNAqc)
library(mobster)


# input parameters
VAF_min <- snakemake@params[["VAF_min"]]
VAF_max <- as.numeric(snakemake@params[["VAF_max"]])
print(paste0("VAF limits: ", VAF_min, ", ", VAF_max))

mobster_cutoff <- as.numeric(snakemake@params[["mobster_cutoff"]])
print(paste0("mobster cutoff: ", mobster_cutoff))

mutations_file = snakemake@input[["SNV_annotated"]] # input SNV data
SNV_ann <- read.table(mutations_file, header = TRUE, fill = TRUE)

print(paste0("Number input SNVs: ", nrow(SNV_ann))) 
print(paste0("Number of NA VAFs: ", sum(is.na(SNV_ann$VAF)))) # dummy check
SNV_ann <- SNV_ann %>% filter(!(is.na(VAF)))

CNA_file <- snakemake@input[["CNA"]] # input CNA data
CNA <- read.table(CNA_file, header = TRUE)
CNA <- CNA %>% select(chromosome, start.pos, end.pos, A, B)
colnames(CNA) <- c("chr", "from", "to", "Major", "minor")
print("CNA data head")
head(CNA, n = 3)

tp_pl_file <- snakemake@input[["purity_and_ploidy"]] # Sequenza purity and ploidy estimates
tp_and_ploidy_sequenza <- read.table(tp_pl_file, header = TRUE)
tp = tp_and_ploidy_sequenza$cellularity[2]
ploidy = tp_and_ploidy_sequenza$ploidy.estimate[2]
print(paste0("Tumor purity estimate: ", tp, "; Tumor ploidy estimate: ", ploidy))

# initiate CNAqc object
x = init(
  mutations = SNV_ann, 
  cna = CNA,
  purity = tp,
  ref = 'hg38')
print(x)

# plot segments data
segments_plot <- plot_segments(x, highlight = c("1:1", "2:1", "2:0", "2:2")) 
ggsave(snakemake@output[["segments"]], segments_plot, units = "cm", create.dir = TRUE,  width = 30, height = 7, dpi = 150)

segments_genome_coverage_perc <- plot_karyotypes(x)
ggsave(snakemake@output[["segments_genome_cov"]], segments_genome_coverage_perc, units = "cm", create.dir = TRUE,  width = 7, height = 7, dpi = 150)

data_hist <- ggpubr::ggarrange(
  plot_data_histogram(x, which = 'VAF'),
  plot_data_histogram(x, which = 'DP'),
  plot_data_histogram(x, which = 'NV'),
  ncol = 3, nrow = 1)
ggsave(snakemake@output[["data_histograms"]], data_hist, units = "cm", create.dir = TRUE,  width = 25, height = 9, dpi = 150)

# peaks analysis
set.seed(0) 
x = analyze_peaks(x)
print(x)

# save peaks analysis results plot
options(repr.plot.width=20, repr.plot.height=5)
peaks_plot <- plot_peaks_analysis(x) + ggtitle(paste0("QC status: ", x$peaks_analysis$QC))
ggsave(snakemake@output[["peaks_QC"]], peaks_plot, units = "cm", create.dir = TRUE,  width = 30, height = 7, dpi = 150)

# select QC passed, 1:1 variants
mut_11_PASS <- x$mutations %>% filter(QC_PASS == TRUE, karyotype == "1:1", VAF < 1) 
print(paste0("Number of mutations after QC pass, karyotype 1:1 and VAF < 1 filtering: ", nrow(mut_11_PASS)))

# select lower VAF cutoff based on input: 
if (VAF_min == "lower_peak"){ # if "lower_peak" find lower peak (set to fixed value of 0.1 if lower peak above 0.25),
    peak <- density(mut_11_PASS$VAF)$x[which.max(density(mut_11_PASS$VAF)$y)]
    if (peak > 0.25){
      print("Single peak VAF distribution expected, setting peak to 0.1")
      peak = 0.1
    }
} else { #  otherwise fixed input cutoff
    peak <- as.numeric(VAF_min)
}

# plot variants with filtering cutoffs and save plot
colors <- c("lower VAF cutoff" = "indianred", "upper VAF cutoff" = "steelblue")
labels1 <- c("lower VAF cutoff" = round(peak, 3), "upper VAF cutoff" = VAF_max)
p1 <- ggplot(mut_11_PASS) + 
    geom_histogram(aes(x = VAF), fill = "grey", color = "darkgrey", binwidth = 0.01) + 
    geom_vline(aes(xintercept = peak, color = "lower VAF cutoff"), linetype = 2, linewidth = 1) + 
    geom_vline(aes(xintercept = VAF_max, color = "upper VAF cutoff"), linetype = 2, linewidth = 1) + 
    theme_bw() +  
    scale_color_manual(name = "VAF cutoffs", values = colors, labels = labels1) + 
    ggtitle(paste0("1:1 karyotype, QC_PASS SNVs\n(n = ", nrow(mut_11_PASS), ")")) + 
    theme(text = element_text(size = 12), 
          plot.title = element_text(size = 12))
ggsave(snakemake@output[["SNVs_QCPASS_11_VAF"]], p1, units = "cm", create.dir = TRUE, width = 12, height = 7, dpi = 150)

# remove low and high VAF variants based on cutoffs
print(paste0("Number of clonal (1:1 genotype) mutations that passed QC and VAF < 1: ", nrow(mut_11_PASS)))
mut_11_PASS_VAFfilt <- mut_11_PASS %>% filter(between(VAF, peak, VAF_max))
print(paste0("Number of clonal (1:1 genotype) mutations that passed QC and VAF filtered: ", nrow(mut_11_PASS_VAFfilt)))

# fit MOBSTER
m_fit = mobster_fit(mut_11_PASS_VAFfilt, 
                    fit.type = "MM",   
                    K = 1:3, 
                    samples = 5, 
                    init = "peaks",
                    tail = c(TRUE, FALSE),
                    epsilon = 1e-10, 
                    maxIter = 300, 
                    seed = 12345, 
                    model.selection = "reICL", 
                    parallel = FALSE, 
                    pi_cutoff = 0.02, 
                    N_cutoff = 10, 
                    description = "MOBSTER")

# print MOBSTER best fit and save
print(m_fit$best)
mobster_best_fit_plot <- plot(m_fit$best, cutoff_assignment = mobster_cutoff)
ggsave(snakemake@output[["mobster_best_fit"]], mobster_best_fit_plot, create.dir = TRUE, device='png', width = 9, height = 7, dpi=150)

# plot and save 5 best MOBSTER fits
model_selection_mobster <- plot_model_selection(m_fit, TOP = 5)
ggsave(snakemake@output[["mobster_best_5_plot"]], model_selection_mobster, create.dir = TRUE, device='png', width = 20, height = 10, dpi=150)

# find clusters based on input cluster assignment cutoff
mobster_clusters <- Clusters(m_fit$best, cutoff_assignment = mobster_cutoff)
print(paste0("Number of C1 cluster mutations: ", nrow(mobster_clusters %>% filter(cluster == "C1"))))

# Choose C1 (clonal) cluster mutations
SNV_tracking <- mobster_clusters |> 
    filter(cluster == "C1") |>
    select(chr, from, to, ref, alt, DP, NV, VAF, is_driver, driver_label, segment_id)
print(head(SNV_tracking, n = 3))
print(paste0("Number of mutations in the final set: ", nrow(SNV_tracking)))

# write results
write.table(SNV_tracking, snakemake@output[["SNVs_for_tracking"]], row.names=FALSE, col.names=TRUE, sep="\t", quote = FALSE)

sink()