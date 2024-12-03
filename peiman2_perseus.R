library(PEIMAN2)
library(tidyverse)

# All plugins can be split into 3 parts
# 1. Reading the command line arguments provided by Perseus and parsing the data.
# 2. Perform the desired functionality on the data.
# 3. Write the results to the expected locations in the Perseus formats.

# 1. Parse command line arguments passed in from Perseus,
# including input file and output file paths.
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 2) {
  stop("Do not provide additional arguments!", call.=FALSE)
}
inFile <- args[1]
outFile <- args[2]


# Use PerseusR to read and write the data in Perseus text format.
library(PerseusR)
mdata <- read.perseus(inFile)

# The mdata object can be easily deconstructed into a number of different
# data frames. Check reference manual or help() for full list.
annotMatrix <- annotCols(mdata)

# 2. Run any kind of analysis on the extracted data.
#library(PEIMAN2)
#prot_list <- unname(sapply(mainMatrix$Uniprot, function(x) strsplit(x, ";")[[1]][1]))
#prot_list <- as.data.frame(prot_list)
#enrich1 <- runEnrichment(protein = prot_list, os.name = 'Homo sapiens (Human)')
sgnf_col <- str_which(colnames(annotMatrix), '_Significant')

df_de <- annotMatrix %>%
  filter(annotMatrix[[sgnf_col]] %in% "+")

prot_ids <- str_extract(df_de$Majority.protein.IDs, '^[^;]*')

results <- runEnrichment(protein = prot_ids, os.name = 'Homo sapiens (Human)')


# 3. Create a matrixData object which can be conveniently written to file
# in the Perseus txt format.
outMdata <- matrixData(main = results[, 2:7], annotCols = results[, c(1, 8)])
write.perseus(outMdata, outFile)
