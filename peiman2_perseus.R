install.packages(c('devtools', 'tidyverse'))

devtools::install_github("pnickchi/PEIMAN2")
devtools::install_github("cox-labs/PerseusR")


library(PEIMAN2)
library(PerseusR)
library(tidyverse)

args = commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  stop("Should provide two arguments: paramFile inFile outFile", call. = FALSE)
}

paramFile <- args[1]
inFile <- args[2]
outFile <- args[3]

parameters <- parseParameters(paramFile)
prot_col <- singleChoiceParamValue(parameters = parameters, name = 'Protein IDs column')
prot_col_psea <- singleChoiceParamValue(parameters = parameters, name = 'Protein IDs column PSEA')
sig_col <- singleChoiceParamValue(parameters = parameters, name = 'Significance column')
sig_col_psea <- singleChoiceParamValue(parameters = parameters, name = 'Significance column PSEA')
scores_col_psea <- singleChoiceParamValue(parameters = parameters, name = 'Protein scores column')
program <- singleChoiceParamValue(parameters = parameters, name = 'Program')
organism <- singleChoiceParamValue(parameters = parameters, name = 'Organism')
organism_psea <- singleChoiceParamValue(parameters = parameters, name = 'Organism PSEA')
adj_method <- singleChoiceParamValue(parameters = parameters, name = 'Adjustment method')
adj_method_psea <- singleChoiceParamValue(parameters = parameters, name = 'Adjustment method PSEA')
pexponent <- intParamValue(parameters = parameters, name = 'Enrichment weighting exponent, p')
nperm <- intParamValue(parameters = parameters, name = 'Number of permutation')
sig_level <- intParamValue(parameters = parameters, name = 'Significance level')
min_size <- intParamValue(parameters = parameters, name = 'Min size')

mdata <- read.perseus(inFile)

mainMatrix <- main(mdata)
annotMatrix <- annotCols(mdata)

if (program == "Singular Enrichment Analysis (SEA)") {
  
  df_de <- annotMatrix %>%
    dplyr::filter(get(sig_col) %in% "+")
  
  prot_ids <- str_extract(df_de[[prot_col]], '^[^;]*')
  
  results <- runEnrichment(protein = prot_ids, os.name = organism, p.adj.method = adj_method)
  
  results <- data.frame(results)
  n_df <- data.frame(n = 1:nrow(results))
  
  outdata <- matrixData(main = n_df, annotCols = results)
  
} else {
  
  sgnf_prots_rows <- which(annotMatrix[[sig_col_psea]] %in% '+')
  
  prot_ids <- str_extract(annotMatrix[[prot_col_psea]][sgnf_prots_rows], '^[^;]*')
  scores <- annotMatrix[[scores_col_psea]][sgnf_prots_rows]
  
  psea_data <- data.frame(UniProtAC = prot_ids,
                          Score = scores)
  
  psea_res <- runPSEA(protein = psea_data,
                      os.name = organism_psea,
                      nperm = nperm,
                      pexponent = pexponent,
                      p.adj.method = adj_method_psea,
                      sig.level = sig_level,
                      minSize = min_size)
  
  res <- data.frame(psea_res[[1]])
  n_df <- data.frame(n = 1:nrow(res))
    
  outdata <- matrixData(main = n_df, annotCols = res)
  
}

write.perseus(outdata, outFile)
