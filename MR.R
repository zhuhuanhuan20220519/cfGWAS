library(TwoSampleMR)
library(dplyr)
library(ieugwasr)

output_file <- "/local_ld_clump/MR_result_local_0.01.csv"
if (!file.exists(output_file)) {
  write.table(
    x = data.frame(), 
    file = output_file,
    sep = ",",
    col.names = TRUE, 
    row.names = FALSE 
  )
}


exposure_dir <- "/cfDNA_GWAS/wuhan_motif_MR/motif_sig"
outcome_dir <- "/cfDNA_GWAS/wuhan_motif_MR/wuhan_phenotype_nosig"

 
exposure_files <- list.files(exposure_dir, pattern = "\\.sig$", full.names = TRUE)

      
outcome_files <- list.files(outcome_dir, pattern = "\\.nosig$", full.names = TRUE)


mr_results_list <- list()

result_list <- list()
r<-data.frame()
# Loop over each exposure file
for (ef in exposure_files) {
  # Read exposure data
  exp_data <- read_exposure_data(ef, sep = "\t")
  exp_data["exposure"]<-sub(".sig","",basename(ef))
  #exp_data <- exp_data %>% filter(!is.na(pval.exposure) & pval.exposure != "")
  if (nrow(exp_data) <= 2) {
    next # Skip to the next file
  }
  clump_data<- ld_clump(dat=dplyr::tibble(rsid=exp_data$SNP,pval=exp_data$pval.exposure,id=exp_data$exposure),
                        clump_kb = 10000, clump_r2 = 0.01, clump_p = 1,
                        bfile="/Downloads/1kg.v3/EAS",
                        plink_bin="/Downloads/plink_mac_20231211/plink")
  if (nrow(clump_data) <= 2) {
    next # Skip to the next file
  }
  #exp_data <- clump_data(exp_data, clump_kb = 10000, clump_r2 = 0.5, clump_p1 = 1, clump_p2 = 1,pop="EAS")
  exp_data<-exp_data %>% filter(SNP %in% clump_data$rsid)
  
  # Loop over each outcome file
  for (of in outcome_files) {
    # Read outcome data
    out_data <- read_outcome_data(of, sep = "\t")
    out_data["outcome"] <- sub(".nosig", "", basename(of))
    
    # Harmonise data
    h <- harmonise_data(exposure_dat = exp_data, outcome_dat = out_data)
    if (nrow(h) <= 2) {
      next # Skip to the next file
    }
    # Check if effect alleles match and adjust if necessary
    h <- h %>%
      mutate(
        beta.exposure = if_else(effect_allele.exposure != effect_allele.outcome, -beta.exposure, beta.exposure),
        effect_allele.exposure = if_else(effect_allele.exposure != effect_allele.outcome, other_allele.exposure, effect_allele.exposure),
        other_allele.exposure = if_else(effect_allele.exposure != effect_allele.outcome, effect_allele.exposure, other_allele.exposure)
      )
    
    # Perform MR analysis
    dat <- mr(h, parameters = default_parameters(), method_list = c("mr_egger_regression", "mr_ivw", "mr_weighted_median"))
    
    basename_prefix<-paste0(sub(".sig", "", basename(ef)),"_",sub(".nosig", "", basename(of)))
    
    exposure <- dat$exposure
    outcome <- dat$outcome
    method <- dat$method
    nsnp <- dat$nsnp
    b <- dat$b
    se <- dat$se
    pval <- dat$pval
    result_df <- data.frame(exposure, outcome, method, nsnp, b, se, pval)
    r<-rbind(r,result_df)
    result_list[[basename_prefix]] <- result_df

    write.table(
      x = result_df,
      file = output_file,
      sep = ",", 
      col.names = FALSE, 
      row.names = FALSE, 
      append = TRUE 
    )
  }
}