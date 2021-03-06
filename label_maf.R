library(GetoptLong)
library(dplyr)
library(tidyr)
library(readr)

## A script to add Tumor/Normal sample barcodes and Study
## to a MAF file.

## Future feature: arbitrary field setting with a -C (custom) flag.

read_maf <- function(fi){
    return(readr::read_tsv(fi,
                           comment="#",
                           col_types=cols(Chromosome=col_character(),
                                          HGNC_Previous_Name=col_character()),
                           skip_empty_rows=TRUE,
                           trim_ws=TRUE
                           ))
}

label_maf <- function(x, var_name, var_value){
    x <- x %>% mutate(!!var_name := var_value)
    return(x)
}

GetoptLong(
    "maf=s", "The MAF file to label.",
    "tumor=s", "The Tumor_Sample_Barcode to add.",
    "normal=s", "The Normal_Sample_Barcode to add.",
    "study=s", "The study field value to add.",
    "caller=s", "The variant caller used to generator the MAF",
    "output=s", "An output file name to write the labeled maf to."
    #"CUSTOM=s%", "A comma-separated list of custom-key-values to add."
)

maf <- read_maf(maf)
maf <- label_maf(maf, "Tumor_Sample_Barcode", tumor)
maf <- label_maf(maf, "Matched_Norm_Sample_Barcode", normal)
maf <- label_maf(maf, "Study", study)
maf <- label_maf(maf, "caller", caller)

#cat(readr::format_tsv(maf))
write_tsv(maf, output)
