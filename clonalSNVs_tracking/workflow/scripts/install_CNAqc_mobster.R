sink(snakemake@log[[1]], append=TRUE) # logging

# Packages
library(devtools)

# install dependencies and packages of interest
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repos = "http://cran.us.r-project.org", INSTALL_opts = '--no-lock')
print("BiocManager installed")
BiocManager::install("BSgenome")
print("BSgenome installed")
devtools::install_github("caravagnalab/CNAqc")
print("CNAqc installed")
devtools::install_github("caravagnalab/mobster")
print("MOBSTER installed")

# confirm installation
f_out <- file(snakemake@output[["installation_check"]])
writeLines(c("Done"), f_out)
close(f_out)

sink()