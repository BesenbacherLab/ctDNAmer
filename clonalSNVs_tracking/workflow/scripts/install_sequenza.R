sink(snakemake@log[[1]], append=TRUE) # logging

# packages
library(remotes)
install_github("aroneklund/copynumber")
remotes::install_bitbucket("sequenzatools/sequenza", dependencies = FALSE)

# confirm installation
f_out <- file(snakemake@output[["installation_check"]])
writeLines(c("Done"), f_out)
close(f_out)

sink()