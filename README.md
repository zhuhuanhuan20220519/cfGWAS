# cfGWAS reveal genetic basis of cell-free DNA features

## Overview
This repository serves as a dedicated space for housing the codebase used in our publication for cfGWAS reveal genetic basis of cell-free DNA features. It is intended to provide researchers with access to the methodologies and algorithms employed in our study, facilitating further research and analysis in the field of genetic.

## License
The code within this repository is licensed under the [MIT License](./LICENSE). Please refer to the license file for more information on the terms and conditions of using and contributing to this project.

## 1. Sequecing data preprocessing and genotype imputaion
### 1.1. Alignment
code：workflow_bwa.sh
```bash
### static path for tools and reference on local disk.
### Tools
gatk=./software/gatk-4.0.4.0/gatk
java=./software/jdk1.8.0_131/bin/java
picard=./software/picard-2.10.10/picard.jar
bwa=./software/bwa/0.7.16/bwa
samtools=./software/samtools
fastp=./software/fastp-0.23.2/fastp

### Reference Genome
hg38=./toolset/reference/hg38_test/Homo_sapiens_assembly38.fasta

### GATK bundle
gatk_bundle_dir=./database/ftp.broadinstitute.org/gsapubftp-anonymous/bundle/hg38
known_indels=./pub/database/ftp.broadinstitute.org/gsapubftp-anonymous/bundle/hg38/beta/Homo_sapiens_assembly38.known_indels.vcf.gz

### Input

sample_id=$1
lane_id=$2
fq=$3
out_path=$4

### Output
final_outdir=$out_path/final_data/$sample_id
outdir=$out_path/temp_data/$sample_id

if [ ! -d $final_outdir ]
then mkdir -p $final_outdir
fi

if [ ! -d $outdir ]
then mkdir -p $outdir
fi

######################################################################################
################################### Pipeline #########################################
######################################################################################
echo "We're doing the job in $sample_id"
echo "We are calculating $fq"
echo "We are doing $lane_id"
echo "We'll save it in $final_outdir"

### step 0: fastp for QC

time $fastp -i $fq -o $outdir/${lane_id}.clean.fq.gz --qualified_quality_phred=5 --unqualified_percent_limit=50 --n_base_limit=10 \
--adapter_sequence="AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA" --adapter_sequence_r2="AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG" \
--disable_trim_poly_g --thread=16 -j $outdir/${lane_id}.json -h $outdir/${lane_id}.html -R $lane_id

### step 1: bwa

echo ""
time $bwa aln -e 10 -t 4 -i 5 -q 0 $hg38 $outdir/${lane_id}.clean.fq.gz > $outdir/${lane_id}.sai && \
    time $bwa samse -r "@RG\tID:${sample_id}\tPL:COMPLETE\tSM:${lane_id}" $hg38 $outdir/${lane_id}.sai $outdir/${lane_id}.clean.fq.gz | $samtools view -h -Sb - > $outdir/${lane_id}.bam && echo "** bwa done **" && \
    time $samtools sort -@ 8 -O bam -o $outdir/${lane_id}.sorted.bam $outdir/${lane_id}.bam && echo "** bam sorted done **" && \
    time $samtools rmdup $outdir/${lane_id}.sorted.bam $outdir/${lane_id}.sorted.rmdup.bam && echo "** rmdup done **" && \
    time $samtools index $outdir/${lane_id}.sorted.rmdup.bam && echo "** index done **" && touch ${outdir}/bwa_sort_rmdup.finish

if [ ! -f ${outdir}/bwa_sort_rmdup.finish ]
then echo "** [WORKFLOW_ERROR_INFO] bwa_sort_rmdup not done **" && exit
fi

### step 2: gatk

echo ""
time $gatk BaseRecalibrator \
    -R $hg38 \
    -I $outdir/${lane_id}.sorted.rmdup.bam \
    --known-sites $gatk_bundle_dir/dbsnp_146.hg38.vcf.gz \
    --known-sites $gatk_bundle_dir/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
    --known-sites $known_indels \
    -O $outdir/${lane_id}.recal_data.table && echo "** BaseRecalibrator done " && touch ${outdir}/baseRecalibrator.finish

if [ ! -f ${outdir}/baseRecalibrator.finish ]
then echo "** [WORKFLOW_ERROR_INFO] baseRecalibrator not done **" && exit
fi

time $gatk ApplyBQSR \
    -R $hg38 \
    --bqsr-recal-file $outdir/${lane_id}.recal_data.table \
    -I $outdir/${lane_id}.sorted.rmdup.bam \
    -O $outdir/${lane_id}.sorted.rmdup.BQSR.bam && echo "** PrintReads done **" && touch ${outdir}/PrintReads.finish

if [ ! -f ${outdir}/PrintReads.finish ]
then echo "** [WORKFLOW_ERROR_INFO] PrintReads not done **" && exit
fi

### step 3: bam index

time $samtools index $outdir/${lane_id}.sorted.rmdup.BQSR.bam && echo "** bam index done **" && touch ${outdir}/bam_index.finish

if [ ! -f ${outdir}/bam_index.finish ]
then echo "** [WORKFLOW_ERROR_INFO] bam index not done **" && exit
fi

### step 4: bam stats
time $samtools stats $outdir/${lane_id}.sorted.rmdup.BQSR.bam > $outdir/${lane_id}.sorted.rmdup.BQSR.bamstats && echo "** bamstats done **" && touch ${outdir}/bamstats.finish

if [ ! -f ${outdir}/bamstats.finish ]
then echo "** [WORKFLOW_ERROR_INFO] bamstats not done **" && exit
fi

### move to final_dir

mv -f $outdir/${lane_id}.sorted.rmdup.BQSR.bam* $final_outdir && echo "** move2final done **" && touch ${outdir}/move2final.finish
if [ ! -f ${outdir}/move2final.finish ]
then echo "** [WORKFLOW_ERROR_INFO] move2final not done **" && exit #v1.6
fi

### clear up
rm -vrf $outdir
```

usage, for example:
``` bash
bash workflow_bwa.sh CL100045254_L01_23 CL100045254_L01_23 CL100045254_L01_23.fq.gz output_path/
```

### 1.2. genotype imputation
code: STITCH.R
```R
#! R
args=commandArgs(T)
outputdir=args[1]
bamlist=args[2]
ref=args[3]
human_posfile=args[4]
human_K=as.numeric(args[5])
human_nGen=as.numeric(args[6])
nCores=as.numeric(args[7])
niterations=as.numeric(args[8])
regionStart=as.numeric(args[9])
regionEnd=as.numeric(args[10])
chr=as.character(args[11])
buffer=as.numeric(args[12])
human_reference_sample_file=as.character(args[13])
human_reference_legend_file=as.character(args[14])
human_reference_haplotype_file=as.character(args[15])
sampleName=as.character(args[16])

tempdir=outputdir
setwd(outputdir)
library("STITCH",lib.loc="./software/Rpackages-3.5.1/")
STITCH(
  bamlist = bamlist,
  reference = ref,
  outputdir = outputdir,
  method = "diploid",
  regenerateInput = TRUE,
  regionStart = regionStart,
  regionEnd = regionEnd,
  buffer = buffer,
  niterations = niterations,
  chr = chr,
  sampleNames_file = sampleName,
  inputBundleBlockSize = 100,
  reference_populations = c("CHB", "CHS", "CDX"),
  reference_haplotype_file = human_reference_haplotype_file,
  reference_sample_file = human_reference_sample_file,
  reference_legend_file = human_reference_legend_file,
  posfile = human_posfile, K = human_K, tempdir = tempdir, nCores = nCores, nGen = human_nGen)

```
usage: run_STITCH.sh
```bash

mkdir ./chr1_1_5000000

time ./software/01.Usr/bin/Rscript \
toolset/GWAS_creatCMD/STITCH.R \
./chr1_1_5000000 \
./bam_forge.list \
./toolset/reference/hg38/Homo_sapiens_assembly38.fasta \
./toolset/reference/1kg.easaf0.01/chr1.pos.txt \
10 16000 4 40 \
1 5000000 chr1 \
250000 \
./toolset/reference/liftover_easaf0.01/1000GP_Phase3.sample \
./toolset/reference/liftover_easaf0.01/chr1.legend.gz \
./toolset/reference/liftover_easaf0.01/chr1.hap.gz \
sampleid_forge.list

rm -rv ./chr1_1_5000000/input \
&& touch ./chr1_1_5000000/input.removed
```

### 1.3. Variants detection
usage: run_basevar.sh
```bash
BaseVarC basetype \
--input bam_forge.list \
--reference Homo_sapiens_assembly38.fasta \
--region chr12:1000001-1500000 \
--output ./chr12-1000001-1500000.10 \
--batch 10 --rerun 
```

### 1.4 Motifs Qunatification

```bash
# code here 
```

## 2. Statistical analyse

### 2.1. Load in Plink & PCA

usage: run_plink_load.sh

```bash
plink2
--const-fid 0 \
--freq \
--hardy \
--make-pgen vzs \
--out WholeGenome2 \
--pca 5 \
--set-all-var-ids chr@:# \
--vcf WG.vcf dosage=DS
```

### 2.2. Genome wide association study (GWAS)
usage: run_GWAS.sh
```bash
./software/plink2 \
--pfile ./plink-format/STITCH/WholeGenome2 'vzs' \
--read-freq ./plink-format/STITCH/WholeGenome2.afreq \
--glm \
--maf 0.05 \
--hwe 1e-6 \
--geno 0.1 dosage \
--pheno ./phenotype_genaral/pheno.csv \
--pheno-quantile-normalize \
--covar ./plink-format/STITCH/basevar_PCA_age.eigenvec \
--covar-variance-standardize \
--out ./output/GENERAL/GWAS
```

### 2.3. Heritability and genetic correlation

code: h2.sh
```bash
glm_add=$1
name=$2
ofile=$3

# mkdir
mkdir ${ofile}/${name}
ofiles=${ofile}/${name}

# glm 2 ss
/share/app/python/3.8.6/bin/python ./toolset/glm_2_ss.py $1 ${ofiles}/${name}.ss

# ss 2 ldsc
/share/app/python/3.8.6/bin/python ./toolset/ss_2_ldsc.py ${ofiles}/${name}.ss ${ofiles}/${name}.ldsc

# get_size
sample_size=$(awk 'NR==2{print $2}' ${ofiles}/${name}.ldsc)
echo "${name}, sample_size: ${sample_size}"

# ldsc 2 munge
./envs/python27/bin/python2.7 \
./software/ldsc-master/munge_sumstats.py \
--sumstats ${ofiles}/${name}.ldsc \
--N $sample_size \
--out ${ofiles}/${name} \
--merge-alleles ./toolset/reference/wuhan_ss.snplist

# munge 2 h2
./envs/python27/bin/python2.7 \
./software/ldsc-master/ldsc.py \
--h2 ${ofiles}/${name}.sumstats.gz \
--ref-ld-chr ./toolset/LDSC/eas_ldscores/ \
--w-ld-chr ./toolset/LDSC/eas_ldscores/ \
--out ${ofiles}/${name}.h2
```
code: rg.sh
```bash
ss1=$1
ss2=$2
ofile=$3

./envs/python27/bin/python2.7 \
./software/ldsc-master/ldsc.py \
--rg $ss1,$ss2 \
--ref-ld-chr ./LDSC/eas_ldscores/ \
--w-ld-chr ./LDSC/eas_ldscores/ \
--out $ofile
```
usage:
```bash
bash glm_2_h2.sh ./output/GENERAL/GWAS.phenotype.add phenotype ./output/GENERAL/h2
bash rg.sh ./output/GENERAL/munge/phenotype1.ldsc.sumstats.gz ./output/GENERAL/munge/phenotype2.ldsc.sumstats.gz ./output/GENERAL/rg/RG_phenotype1_phenotype2
```

### 2.4. **Heritability partition**

code: ph2.sh

```bash
./envs/python27/bin/python2.7 \
./software/ldsc-master/ldsc.py \
--h2-cts ${ofiles}/${name}.sumstats.gz \
--ref-ld-chr ./ldsc_reference/1000G_Phase3_EAS_baselineLD_v2.2_ldscore/baselineLD. \
--w-ld-chr ./ldsc_reference/1000G_Phase3_EAS_weights_hm3_no_MHC/weights.EAS.hm3_noMHC. \
--frqfile-chr ./ldsc_reference/1000G_Phase3_EAS_plinkfiles/1000G.EAS.QC. \
--ref-ld-chr-cts ./ldsc_reference/Multi_tissue_gene_expr.EAS.ldcts \
--out ${ofiles}/${name}.ph2
```

### 2.5. **Pathway-based analysis**

code: PASCAL.sh

```bash
cd ./PASCAL/
./PASCAL/Pascal \
--pval=${ifile} \
--runpathway=on \
--custom=ASN \
--customdir=./PASCAL/resources/ASN/CHR \
--outsuffix=${ofile}
```


### 2.6. Mendelian randomization analysis
code: MR.R
```R
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

```



### 2.7. Colocalization analysis

code: coloc.R

```R

library("coloc")
library("ggplot2")
library("locuscomparer")
library(data.table)

setwd("D:\\cfDNA\\coloc\\")
rs <- fread("chrall.arrange.1000genomes.chrpos_2_rs.index")
rs=rs[,3:4]
colnames(rs)[2] <- "ID"

wh_folder <- "D:\\cfDNA\\coloc\\new_coloc\\plot_data\\WH"  # Data 1 folder path
motif_folder <-  "D:\\cfDNA\\coloc\\new_coloc\\plot_data\\motif"  # Data 2 folder path

wh_files <- list.files(wh_folder, pattern = '*.txt', full.names = TRUE)  # Get all files in the data1 folder
motif_files <- list.files(motif_folder,pattern = '*.txt', full.names = TRUE)  # Get all files in the data2 folder

output <- data.frame(gene = character(),
                     trait1 = character(),
                     trait2 = character(),
                     PP0  = numeric(),
                     PP1  = numeric(),
                     PP2  = numeric(),
                     PP3  = numeric(),
                     PP4  = numeric(),
                     nsnps = numeric()
                     )

colnames(output) <- c("gene", "trait1", "trait2","PP0", "PP1", "PP2","PP3","PP4","nsnps")
##coloc
for (i in 1:length(wh_files)) {
  wh_prefix <- tools::file_path_sans_ext(basename(wh_files[i]))  
  gene_prefix <- gsub("_.*","", wh_prefix) 
  motif_matches <- grep(gene_prefix, motif_files, value = TRUE)  
  ##Read file 1
  wh <- read.table(file=wh_files[i], header=F, as.is=T)
  colnames(wh) <- c("trait","CHROM","POS","SNP","GENE","REGION","A2","A1","MAF","OBS_CT","BETA","SE","zsores","P")
  wh$ID <- paste(wh$ID, wh$CHROM, wh$POS, sep=";")
  ## Check if there is a row in column P that is less than 5e-8. If not, exit the loop.
  if(all(wh$P >= 5e-8)){next
  }
  for (j in 1:length(motif_matches)) {
    ##Read file 2
    motif <- read.table(file=motif_matches[j], header=F, as.is=T)
    colnames(motif) <- c("CHROM","POS","ID","REF","ALT","A1","TEST","OBS_CT","BETA","SE","T_STAT","P")
    ##Merge file 1 file 2
    input <- merge(wh,motif, by="ID", all=FALSE, suffixes=c("_wh","_motif"))
    ##Corrected beta value
    input$BETAN_wh <- ifelse(input$A1_wh==input$A1_motif, input$BETA_wh, -input$BETA_wh)
    ##Calculate variance
    input$VAR_wh <- input$SE_wh*input$SE_wh
    input$VAR_motif <- input$SE_motif*input$SE_motif
    
    # Perform coloc analysis
    result <- coloc.abf(dataset1=list(snp=input$ID,pvalues=input$P_wh, beta=input$BETAN_wh,varbeta=input$VAR_wh,type="quant", N=input$OBS_CT_wh[1]), 
                        dataset2=list(snp=input$ID,pvalues=input$P_motif,beta=input$BETA_motif,varbeta=input$VAR_motif,type="quant", N=input$OBS_CT_motif[1]),
                        MAF = input$MAF,p12 = 1e-5)
    
    # Save the results
    gene <- gene_prefix
    trait1 <- gsub(".*_","", wh_prefix)
    trait2 <- gsub(".*_([^.]*)\\..*", "\\1", tools::file_path_sans_ext(basename(motif_matches[j])))
    PP0 <- result[["summary"]][["PP.H0.abf"]]
    PP1 <- result[["summary"]][["PP.H1.abf"]]
    PP2 <- result[["summary"]][["PP.H2.abf"]]
    PP4 <- result[["summary"]][["PP.H4.abf"]]
    PP3 <- result[["summary"]][["PP.H3.abf"]]
    nsnps <- result[["summary"]][["nsnps"]]
    output <- rbind(output,c(gene, trait1, trait2,PP0, PP1, PP2,PP3,PP4,nsnps))
    }
}

write.table(output,"/cfDNA/coloc/output/new_coloc_result1e-5.txt",sep = "\t",quote = "F",row.names = F)



##locuszoom
for (i in 1:length(wh_files)) {
  wh_prefix <- tools::file_path_sans_ext(basename(wh_files[i]))  
  gene_prefix <- gsub("_.*","", wh_prefix) 
  motif_matches <- grep(gene_prefix, motif_files, value = TRUE)  
  for (j in 1:length(motif_matches)) {
    ##Read file 1
    print(i)
    print(j)
    trait1 <- gsub(".*_","", wh_prefix)
    trait2 <- gsub(".*_([^.]*)\\..*", "\\1", tools::file_path_sans_ext(basename(motif_matches[j])))
    wh <- read.table(file=wh_files[i], header=F, as.is=T)
    colnames(wh) <- c("CHROM","POS","ID","REF","ALT","A1","TEST","OBS_CT","BETA","SE","T_STAT","P")
    
    ##Read file 2
    motif <- read.table(file=motif_matches[j], header=F, as.is=T)
    colnames(motif) <- c("CHROM","POS","ID","REF","ALT","A1","TEST","OBS_CT","BETA","SE","T_STAT","P")
    
    input <- merge(wh,motif, by="ID", all=FALSE, suffixes=c("_wh","_motif"))
    
    input2 <- merge(input,rs,by='ID',all.x = TRUE)
    write.table(input2,paste0("plot//",gene_prefix,".txt"),row.names = F,sep = "\t",quote = F)

  }}



library(biomaRt)
snpmart = useMart(biomart = "ENSEMBL_MART_SNP", dataset="hsapiens_snp")
plotdata <- list.files("D:\\cfDNA\\coloc\\new_coloc\\plot_data\\plot\\", pattern = '*.txt', full.names = TRUE)
for (i in 1:length(plotdata)){
  input2 <- read.csv(plotdata[i],sep = "\t")
  df <- input2[is.na(input2$SNP),]

  chr <- df$CHROM_wh
  pos <- df$POS_wh

  info <- data.frame()

  for (j in 1:length(chr)){
    snpinfo <- getBM(attributes = c('refsnp_id','allele','chr_name','chrom_start','chrom_end'), 
                  filters = c('chr_name','start','end'), 
                   values = list(chr[j],pos[j],pos[j]), 
                   mart = snpmart)
    snpinfo <- snpinfo[snpinfo$chrom_start==snpinfo$chrom_end,]
    snpinfo$ID <-  paste0("chr",snpinfo$chr_name,":",snpinfo$chrom_start)
    snpinfo <- snpinfo[,c("ID","refsnp_id","allele")]
    info <- rbind(info, snpinfo)
  }

  info <- info[duplicated(info$ID)==F,]
  input2 <- merge(input2,info,by = "ID",all.x = T)
  input2$SNP2 <- ifelse(is.na(input2$SNP),input2$refsnp_id,input2$SNP)

  gene_prefix <- tools::file_path_sans_ext(basename(plotdata[i]))
  trait1 <- gsub(".*_","", tools::file_path_sans_ext(basename(wh_files[i])))
  trait2 <- gsub(".*_", "", tools::file_path_sans_ext(basename(motif_files[i])))

  l1 <- input2[, c("SNP2", "P_wh")]
  colnames(l1) <- c("rsid","pval")
  l2 <- input2[,c("SNP2","P_motif")]
  colnames(l2) <- c("rsid","pval")

  p1 <- locuscompare(in_fn1=l1, in_fn2=l2, title1=paste0("GWAS_",trait1), title2=paste0("GWAS_",trait2),
                       population = "EAS",
                       genome="hg38",
    )

  print(paste0("D:\\cfDNA\\coloc\\new_coloc\\plot\\",gene_prefix,"_",trait1,"_",trait2))
  ggsave(paste0("D:\\cfDNA\\coloc\\new_coloc\\plot\\",gene_prefix,"_",trait1,"_",trait2,".png"),plot =p1, height = 5, width = 10,dpi=200)
}

```
