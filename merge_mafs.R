library(dplyr)
library(tidyr)
library(tibble)
library(readr)
library(stringr)
library(GetoptLong)

## A script to add Tumor/Normal sample barcodes and Study
## to a MAF file.

## Future feature: arbitrary field setting with a -C (custom) flag.
read_maf <- function(fi){
    return(readr::read_tsv(fi,
                           comment="#",
                           col_types=cols(Chromosome=col_character(),
                                          HGNC_Previous_Name=col_character(),
                                          dbSNP_RS=col_character(),
                                          dbSNP_ASP=col_character(),
                                          dbSNP_ASS=col_character(),
                                          dbSNP_CDA=col_character(),
                                          dbSNP_CFL=col_character(),
                                          dbSNP_DSS=col_character(),
                                          dbSNP_G5=col_character(),
                                          CGC_Mutation_Type=col_character(),
                                          `HGNC_OMIM_ID(supplied_by_OMIM)`=col_character()),
                           skip_empty_rows=TRUE,
                           trim_ws=TRUE
                           ))
}

read_funcotator <- function(fi){
    return(readr::read_tsv(fi,
                           skip_empty_rows=TRUE,
                           trim_ws=TRUE,
                           comment="#",
                           col_types=cols(
                                          Hugo_Symbol=col_character(),
                                          Entrez_Gene_Id=col_character(),
                                          Center=col_character(),
                                          NCBI_Build=col_character(),
                                          Chromosome=col_character(),
                                          Start_Position=col_integer(),
                                          End_Position=col_integer(),
                                          Strand=col_character(),
                                          Variant_Classification=col_character(),
                                          Variant_Type=col_character(),
                                          Reference_Allele=col_character(),
                                          Tumor_Seq_Allele1=col_character(),
                                          Tumor_Seq_Allele2=col_character(),
                                          dbSNP_RS=col_character(),
                                          dbSNP_Val_Status=col_character(),
                                          Tumor_Sample_Barcode=col_character(),
                                          Matched_Norm_Sample_Barcode=col_character(),
                                          Match_Norm_Seq_Allele1=col_character(),
                                          Match_Norm_Seq_Allele2=col_character(),
                                          Tumor_Validation_Allele1=col_character(),
                                          Tumor_Validation_Allele2=col_character(),
                                          Match_Norm_Validation_Allele1=col_character(),
                                          Match_Norm_Validation_Allele2=col_character(),
                                          Verification_Status=col_character(),
                                          Validation_Status=col_character(),
                                          Mutation_Status=col_character(),
                                          Sequencing_Phase=col_character(),
                                          Sequence_Source=col_character(),
                                          Validation_Method=col_character(),
                                          Score=col_integer(),
                                          BAM_File=col_character(),
                                          Sequencer=col_character(),
                                          Tumor_Sample_UUID=col_character(),
                                          Matched_Norm_Sample_UUID=col_character(),
                                          Genome_Change=col_character(),
                                          Annotation_Transcript=col_character(),
                                          Transcript_Strand=col_character(),
                                          Transcript_Exon=col_character(),
                                          Transcript_Position=col_character(),
                                          cDNA_Change=col_character(),
                                          Codon_Change=col_character(),
                                          Protein_Change=col_character(),
                                          Other_Transcripts=col_character(),
                                          Refseq_mRNA_Id=col_character(),
                                          Refseq_prot_Id=col_character(),
                                          SwissProt_acc_Id=col_character(),
                                          SwissProt_entry_Id=col_character(),
                                          Description=col_character(),
                                          UniProt_AApos=col_character(),
                                          UniProt_Region=col_character(),
                                          UniProt_Site=col_character(),
                                          UniProt_Natural_Variations=col_character(),
                                          UniProt_Experimental_Info=col_character(),
                                          GO_Biological_Process=col_character(),
                                          GO_Cellular_Component=col_character(),
                                          GO_Molecular_Function=col_character(),
                                          COSMIC_overlapping_mutations=col_integer(),
                                          COSMIC_fusion_genes=col_character(),
                                          COSMIC_tissue_types_affected=col_character(),
                                          COSMIC_total_alterations_in_gene=col_integer(),
                                          Tumorscape_Amplification_Peaks=col_character(),
                                          Tumorscape_Deletion_Peaks=col_character(),
                                          TCGAscape_Amplification_Peaks=col_character(),
                                          TCGAscape_Deletion_Peaks=col_character(),
                                          DrugBank=col_character(),
                                          ref_context=col_character(),
                                          gc_content=col_double(),
                                          CCLE_ONCOMAP_overlapping_mutations=col_integer(),
                                          CCLE_ONCOMAP_total_mutations_in_gene=col_integer(),
                                          CGC_Mutation_Type=col_character(),
                                          CGC_Translocation_Partner=col_character(),
                                          CGC_Tumor_Types_Somatic=col_character(),
                                          CGC_Tumor_Types_Germline=col_character(),
                                          CGC_Other_Diseases=col_character(),
                                          DNARepairGenes_Activity_linked_to_OMIM=col_character(),
                                          FamilialCancerDatabase_Syndromes=col_character(),
                                          MUTSIG_Published_Results=col_character(),
                                          OREGANNO_ID=col_character(),
                                          OREGANNO_Values=col_character(),
                                          tumor_f=col_double(),
                                          t_alt_count=col_integer(),
                                          t_ref_count=col_integer(),
                                          n_alt_count=col_integer(),
                                          n_ref_count=col_integer(),
                                          Gencode_19_secondaryVariantClassification=col_character(),
                                          Achilles_Top_Genes=col_character(),
                                          CGC_Name=col_character(),
                                          CGC_GeneID=col_character(),
                                          CGC_Chr=col_character(),
                                          CGC_Chr_Band=col_character(),
                                          CGC_Cancer_Somatic_Mut=col_character(),
                                          CGC_Cancer_Germline_Mut=col_character(),
                                          CGC_Cancer_Syndrome=col_character(),
                                          CGC_Tissue_Type=col_character(),
                                          CGC_Cancer_Molecular_Genetics=col_character(),
                                          CGC_Other_Germline_Mut=col_character(),
                                          ClinVar_HGMD_ID=col_character(),
                                          ClinVar_SYM=col_character(),
                                          ClinVar_TYPE=col_character(),
                                          ClinVar_ASSEMBLY=col_character(),
                                          ClinVar_rs=col_character(),
                                          CosmicFusion_fusion_id=col_character(),
                                          DNARepairGenes_Chromosome_location_linked_to_NCBI_MapView=col_character(),
                                          DNARepairGenes_Accession_number_linked_to_NCBI_Entrez=col_character(),
                                          Familial_Cancer_Genes_Synonym=col_character(),
                                          Familial_Cancer_Genes_Reference=col_character(),
                                          Gencode_XHGNC_hgnc_id=col_character(),
                                          HGNC_HGNC_ID=col_character(),
                                          HGNC_Status=col_character(),
                                          HGNC_Locus_Type=col_character(),
                                          HGNC_Locus_Group=col_character(),
                                          HGNC_Previous_Symbols=col_character(),
                                          HGNC_Previous_Name=col_character(),
                                          HGNC_Synonyms=col_character(),
                                          HGNC_Name_Synonyms=col_character(),
                                          HGNC_Chromosome=col_character(),
                                          HGNC_Date_Modified=col_character(),
                                          HGNC_Date_Symbol_Changed=col_character(),
                                          HGNC_Date_Name_Changed=col_character(),
                                          HGNC_Accession_Numbers=col_character(),
                                          HGNC_Enzyme_IDs=col_character(),
                                          HGNC_Ensembl_Gene_ID=col_character(),
                                          HGNC_Pubmed_IDs=col_character(),
                                          HGNC_RefSeq_IDs=col_character(),
                                          HGNC_Gene_Family_ID=col_character(),
                                          HGNC_Gene_Family_Name=col_character(),
                                          HGNC_CCDS_IDs=col_character(),
                                          HGNC_Vega_ID=col_character(),
                                          `HGNC_OMIM_ID(supplied_by_OMIM)`=col_character(),
                                          `HGNC_RefSeq(supplied_by_NCBI)`=col_character(),
                                          `HGNC_UniProt_ID(supplied_by_UniProt)`=col_character(),
                                          `HGNC_Ensembl_ID(supplied_by_Ensembl)`=col_character(),
                                          `HGNC_UCSC_ID(supplied_by_UCSC)`=col_character(),
                                          Oreganno_Build=col_character(),
                                          Simple_Uniprot_alt_uniprot_accessions=col_character(),
                                          chr1_a_bed_name=col_character(),
                                          chr1_a_bed_score=col_character(),
                                          chr1_a_bed_strand=col_character(),
                                          chr1_a_bed_thickStart=col_character(),
                                          chr1_a_bed_thickEnd=col_character(),
                                          chr1_a_bed_itemRgb=col_character(),
                                          chr1_a_bed_blockCount=col_character(),
                                          chr1_a_bed_blockSizes=col_character(),
                                          chr1_a_bed_blockStarts=col_character(),
                                          chr1_b_bed_name=col_character(),
                                          chr1_b_bed_score=col_character(),
                                          chr1_b_bed_strand=col_character(),
                                          chr1_b_bed_thickStart=col_character(),
                                          chr1_b_bed_thickEnd=col_character(),
                                          chr1_b_bed_itemRgb=col_character(),
                                          chr1_b_bed_blockCount=col_character(),
                                          chr1_b_bed_blockSizes=col_character(),
                                          chr1_b_bed_blockStarts=col_character(),
                                          dbSNP_ASP=col_character(),
                                          dbSNP_ASS=col_character(),
                                          dbSNP_CAF=col_character(),
                                          dbSNP_CDA=col_character(),
                                          dbSNP_CFL=col_character(),
                                          dbSNP_COMMON=col_character(),
                                          dbSNP_DSS=col_character(),
                                          dbSNP_G5=col_character(),
                                          dbSNP_G5A=col_character(),
                                          dbSNP_GENEINFO=col_character(),
                                          dbSNP_GNO=col_character(),
                                          dbSNP_HD=col_character(),
                                          dbSNP_INT=col_character(),
                                          dbSNP_KGPhase1=col_character(),
                                          dbSNP_KGPhase3=col_character(),
                                          dbSNP_LSD=col_character(),
                                          dbSNP_MTP=col_character(),
                                          dbSNP_MUT=col_character(),
                                          dbSNP_NOC=col_character(),
                                          dbSNP_NOV=col_character(),
                                          dbSNP_NSF=col_character(),
                                          dbSNP_NSM=col_character(),
                                          dbSNP_NSN=col_character(),
                                          dbSNP_OM=col_character(),
                                          dbSNP_OTH=col_character(),
                                          dbSNP_PM=col_character(),
                                          dbSNP_PMC=col_character(),
                                          dbSNP_R3=col_character(),
                                          dbSNP_R5=col_character(),
                                          dbSNP_REF=col_character(),
                                          dbSNP_RV=col_character(),
                                          dbSNP_S3D=col_character(),
                                          dbSNP_SAO=col_character(),
                                          dbSNP_SLO=col_character(),
                                          dbSNP_SSR=col_character(),
                                          dbSNP_SYN=col_character(),
                                          dbSNP_TOPMED=col_character(),
                                          dbSNP_TPA=col_character(),
                                          dbSNP_U3=col_character(),
                                          dbSNP_U5=col_character(),
                                          dbSNP_VC=col_character(),
                                          dbSNP_VP=col_character(),
                                          dbSNP_WGT=col_character(),
                                          dbSNP_WTD=col_character(),
                                          dbSNP_dbSNPBuildID=col_character(),
                                          dbSNP_ID=col_character(),
                                          dbSNP_FILTER=col_character(),
                                          `HGNC_Entrez_Gene_ID(supplied_by_NCBI)`=col_character(),
                                          dbSNP_RSPOS=col_character(),
                                          dbSNP_VLD=col_character(),
                                          QSI=col_character(),
                                          TQSI=col_character(),
                                          NT=col_character(),
                                          QSI_NT=col_character(),
                                          TQSI_NT=col_character(),
                                          SGT=col_character(),
                                          RU=col_character(),
                                          RC=col_character(),
                                          IC=col_character(),
                                          IHP=col_character(),
                                          MQ=col_integer(),
                                          MQ0=col_integer(),
                                          SOMATIC=col_character(),
                                          OVERLAP=col_character(),
                                          SomaticEVS=col_character(),
                                          Study=col_character(),
                                          caller=col_character()
                                          )
                           ))
}

add_maf_value <- function(x, var_name, var_value){
    x <- x %>% mutate(!!var_name := var_value)
    return(x)
}

maf_cols <- c(
              "Hugo_Symbol",
              "Entrez_Gene_Id",
              "Center",
              "NCBI_Build",
              "Chromosome",
              "Start_Position",
              "End_Position",
              "Strand",
              "Variant_Classification",
              "Variant_Type",
              "Reference_Allele",
              "Tumor_Seq_Allele1",
              "Tumor_Seq_Allele2",
              "dbSNP_RS",
              "dbSNP_Val_Status",
              "Tumor_Sample_Barcode",
              "Matched_Norm_Sample_Barcode",
              "Match_Norm_Seq_Allele1",
              "Match_Norm_Seq_Allele2",
              "Tumor_Validation_Allele1",
              "Tumor_Validation_Allele2",
              "Match_Norm_Validation_Allele1",
              "Match_Norm_Validation_Allele2",
              "Verification_Status",
              "Validation_Status",
              "Mutation_Status",
              "Sequencing_Phase",
              "Sequence_Source",
              "Validation_Method",
              "Score",
              "BAM_File",
              "Sequencer",
              "Tumor_Sample_UUID",
              "Matched_Norm_Sample_UUID"
              )

merge_vars <- c(
                "Chromosome",
                "Start_Position",
                "End_Position",
                "Reference_Allele",
                "Tumor_Seq_Allele1",
                "Tumor_Seq_Allele2",
                "Tumor_Sample_Barcode",
                "Matched_Norm_Sample_Barcode"
                )

relabel_vars <- c(
                  "Chromosome",
                  "Start_Position",
                  "End_Position",
                  "Reference_Allele",
                  "Tumor_Seq_Allele1",
                  "Tumor_Seq_Allele2",
                  "Tumor_Sample_Barcode",
                  "Matched_Norm_Sample_Barcode",
                  "caller"
                  )


merge_maf_dfs <- function(x, z){
    m <- dplyr::full_join(x, z, by = merge_vars)
    return (m)
}

minimize_maf <- function(x){
    mini_maf <- x %>% 
    dplyr::select(Chromosome,
                                    Start_Position,
                                    End_Position,
                                    Reference_Allele,
                                    Tumor_Seq_Allele1,
                                    Tumor_Seq_Allele2,
                                    Tumor_Sample_Barcode,
                                    Matched_Norm_Sample_Barcode,
                                    caller) %>%
    distinct()

    return (mini_maf)
}

convert_caller_column <- function(x){
    return 
    (
     x %>% mutate(ones=1) %>% spread(caller, ones)
     )

}

set_caller_column <- function(x){
    if (! "strelka1" %in% names(x)){
        x <- x %>% mutate(strelka1=0)
    }
    if (! "mutect2" %in% names(x)){
        x <- x %>% mutate(mutect2=0)
    }
    if (! "mutect1" %in% names(x)){
        x <- x %>% mutate(mutect1=0)
    }
    if (! "strelka2" %in% names(x)){
        x <- x %>% mutate(strelka2=0)
    }

    y <- x %>%
    mutate(caller = 
           case_when(
                     strelka1 == 1 ~ "strelka1",
                     strelka2 == 1 ~ "strelka2",
                     mutect1 == 1 ~ "mutect1",
                     mutect2 == 1 ~ "mutect2",
                     TRUE ~ NA_character_
                     )
           )
    return(y)
}

threads=1
output="merged.maf"
m_dframe <- tibble(
                   Chromosome = character(),
                   Start_Position = integer(),
                   End_Position = integer(),
                   Reference_Allele = character(),
                   Tumor_Seq_Allele1 = character(),
                   Tumor_Seq_Allele2 = character(),
                   Tumor_Sample_Barcode = character(),
                   Matched_Norm_Sample_Barcode = character(),
                   caller = character()
                   )

GetoptLong(

           "maf=s@", "MAF file(s). May be repeated multiple times.",
           "output=s","The name of the output file. (Default: merged.maf)",
           "threads=i", "The number of OMP threads for FST compression/decompression.",
           "writefst", "Write an fst version of the final merged MAF."
           )

fst::threads_fst(threads)

for (fi in maf){
    message("Processing: ", fi)
    dname <- dirname(fi)
    bname <- basename(fi)
    m <- read_funcotator(fi)
    fst::write.fst(m, paste(bname, "fst", sep="."), 80)
    mini_maf <- minimize_maf(m)
    mini_maf <- convert_caller_column(mini_maf)
    m_dframe <- merge_maf_dfs(m_dframe, mini_maf)
    rm(m)
}

m_dframe <- set_caller_column(m_dframe)
m_dframe <- m_dframe %>% distinct()

message("cols:", paste(names(m_dframe), sep=" "))

out_dframe <-tibble(
                    Chromosome = character(),
                    Start_Position = integer(),
                    End_Position = integer(),
                    Reference_Allele = character(),
                    Tumor_Seq_Allele1 = character(),
                    Tumor_Seq_Allele2 = character(),
                    Tumor_Sample_Barcode = character(),
                    Matched_Norm_Sample_Barcode = character(),
                    caller = character()
                    )


for (fi in maf){
    bname <- basename(fi)
    fast_maf <- fst::read.fst(paste(bname, "fst", sep=".")) %>% distinct()
    ## Process only one caller's results
    caller_val <- (fast_maf %>% select(caller) %>% distinct())$caller[1]
    caller_dframe <- m_dframe %>% filter(caller == caller_val)
    merged_dframe <- left_join(caller_dframe, fast_maf, by=relabel_vars) 
    out_dframe <- bind_rows(out_dframe, merged_dframe)
}
out_dframe <- out_dframe %>% select(all_of(maf_cols), everything())

if (writefst){
    fst::write.fst(out_dframe, paste(output, "fst", sep="."), 100)
}

write_tsv(out_dframe, output)


