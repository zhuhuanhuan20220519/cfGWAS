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
code: Create_bash_STITCH.py
```python

#!/usr/bin/env python
# coding: utf-8
# lilinxuan@genomics.cn
# stitch for imputation
# v2.0

import os

#Specify the arguments before excute this code
#this is the chunk size, the smaller the size, lower compute resoure required
chunk_size = int('5000000')

#the folder of working
outdir='/path_to/STITCH_analysis'

#sample names and bam paths, need to match each other
sample_namelist='/path_to/STITCH_analysis/bampath_4659.list'
bamlist='/path_to/STITCH_analysis/bamid.list'

#the length of each chromosomes, example file was in GRCh38, please change if running on other versions of reference genome
CHR_length='/path_to/STITCH_analysis/CHR.len'

#the path to the STITCH code
STITCH_script ='/path_to/STITCH_analysis/2.2.run_STITCH.sh'


######################################
###############working################
######################################
#Check if file exist
flag=8*os.path.isfile(bamlist)+4*os.path.isfile(sample_namelist)+2*os.path.isfile(CHR_length)
if flag== 14:
    print("\nUsing file:\nbamlist: "+bamlist +'\nsample_name: '+sample_namelist+'\nCHR.len: '+CHR_length)
else:
    raise ValueError('Some required file (bamlist, sample namelist, etc.) does not exist, please check.')

#mkdir
output_path = outdir + '/output/'

print('\nlog file will write at '+output_path+'/2.1.Create_bash_STITCH.log')

if os.path.isdir(output_path) == False:
    os.mkdir(output_path)
    log = open(output_path+'/2.1.Create_bash_STITCH.log','w')
    log.write('mkdir: '+output_path+'\n')
else:
    log = open(output_path+'/2.1.Create_bash_STITCH.log','w')

log.write('chr_len_path:{0}\noutdir:{1}\nbamlist:{2}\nbin={3}\n'.format(CHR_length,output_path,bamlist,chunk_size))

sh_path = outdir + '/bash/'
if os.path.isdir(sh_path) == False:
    os.mkdir(sh_path)
    log.write('mkdir: '+sh_path+'\n')
else:
    pass

#Create the bash files
chr_len=open(CHR_length ,'r').read().split('\n')[:-1]
chr_len_dic={}
for chrm in chr_len:
    chr_ ,len_=chrm.split(' ')
    chr_len_dic[chr_]=int(len_)

for key,value in chr_len_dic.items():
    if (key == 'chrY')|(key == 'chrM'):
        continue
        #Notice: We ignored chromosome Y and Mitochondrial genome, since the imputation on these genomes are complex
    else:
        sh_this=open(sh_path+key+'.sh','w')

    for i in range(0,1000):#We regulated the number of chunks on a single chromsome with < 1000
        start=1+chunk_size*i
        end=chunk_size*(i+1)
        if(end>value):#end of chr
            sh_this.write(
                f'bash {STITCH_script} {key}_{start}_{end} {key} {start} {value}'
            )
            break
        else:
            sh_this.write(
                f'bash {STITCH_script} {key}_{start}_{end} {key} {start} {end}\n'
            )
    sh_this.close()
log.close()

print('All file prepared.')
```
code: run_STITCH.sh

```bash
#Specify the arguments before excute the code related
#path for bcftools
bcftools=/path_to/bcftools

##sample names and bam paths, the sample file as 2.1.Create_bash_STITCH.py
sample_namelist=/path_to/STITCH_analysis/bampath.list
bamlist=/path_to/STITCH_analysis/bamid.list

#output folder
outdir=/path_to/STITCH_analysis/output

#Reference Genome, should be the same as 01.workflow_bwa.sh
hg38=/path_to/Homo_sapiens_assembly38.fasta

#Before starting STITCH, please make sure you have correctly installed STITCH in your environment
#If correctly installed, you can find the STITCH.R in the directory you install STITCH
STITCH_R=/path_to/STITCH.R

chr_region=$1
chr=$2
start=$3
end=$4

if [ ! -d ${outdir}/${chr}/${chr_region} ]
then mkdir -p ${outdir}/${chr}/${chr_region}
fi

###Specify the arguments before execute!
time /path_to/Rscript $STITCH_R \
--outputdir ${outdir}/${chr}/${chr_region} \
--bamlist $bamlist \
--sampleNames_file $sample_namelist \
--reference $hg38 \
--K 10 --nGen 16000 --nCores 1 \
--regionStart $start --regionEnd $end --chr $chr \
--buffer 250000 
#--niterations 40 \
#--posfile /path_to/${chr}.pos.txt \
#--reference_sample_file /path_to/eas.${chr}.samples \
#--reference_legend_file /path_to/eas.${chr}.legend.gz \
#--reference_haplotype_file /path_to/eas.${chr}.hap.gz


$bcftools index -t ${outdir}/${chr}/${chr_region}/${chr}/${chr_region}/stitch.${chr}.${start}.${end}.vcf.gz

rm -rv  ${outdir}/${chr}/${chr_region}/input \
&& touch  ${outdir}/${chr}/${chr_region}/input.removed
```

code: Merge.sh

```bash
#Specify the arguments before excute the code related
#path for bcftools
bcftools=/path_to/bcftools

$bcftools concat \
-f /path_to/STITCH_analysis/completed_vcf.list \
-Oz \
-o /path_to/STITCH_analysis/final/STITCH.vcf.gz
```

usage, for example:

```bash
python Create_bash_STITCH.py
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

## 3. Motif analyse

code: bam to motif freq.pl

```perl
#!/usr/bin/perl -w
use strict;
use warnings;
use List::Util qw/min max/;

my $samtools="samtools";    #path to samtools
my $fa = "hg38.p13.fa";     # reference genome
my $motif_list = "4mer.motif.list";       #list of 256 types of 4mer motif

# Step 1: Filter raw BAM file to obtain clean SAM file
# Usage: perl script.pl <input.bam> <min_map_qual> <max_mismatch> <motif_length> <sample_prefix>

if (@ARGV != 5) {
    die "Usage: perl script.pl <input.bam> <min_map_qual> <max_mismatch> <motif_length> <sample_prefix>\n";
}

my ($bam_file, $min_map_qual, $max_mismatch, $motif_length, $sample_prefix) = @ARGV;
my %hash;
my %hash_pair;

# Open BAM file and filter reads based on quality, mismatches, and alignment flags
open my $bam_in, "-|", "$samtools view $bam_file" or die "Could not open BAM file: $!\n";

while (<$bam_in>) {
    chomp;
    my $line = $_;
    my @F = split /\s+/, $line;

    # Apply filtering conditions
    next if ($F[6] ne "=");                  # Only keep pairs mapped to the same chromosome
    next if ($F[4] < $min_map_qual);         # Mapping quality cutoff
    next if ($F[1] & 0x0400);                # Exclude PCR/optical duplicates
    next if ($F[1] & 0x4);                   # Exclude unmapped segments
    next if ($F[1] & 0x8);                   # Exclude pairs with unmapped segments
    next if ($F[1] & 0x10 && $F[1] & 0x20);  # Exclude if both reads are reverse complemented
    next unless ($F[1] & 0x10 || $F[1] & 0x20); # Keep if only one read is reverse complemented
    next if ($F[1] & 0x100);                 # Exclude secondary alignments
    next if ($F[8] > 600 || $F[8] < -600);   # Exclude fragments >600 bp in length
    next if (exists $hash{$F[0]});           # Exclude fragments with multiple mismatches

    # Store pairs for further filtering and check mismatches
    $hash_pair{$F[0]}++;
    if ($hash_pair{$F[0]} == 2) {
        delete $hash_pair{$F[0]};
    }

    my $distance = 10;  # Initialize mismatch distance
    if ($F[5] !~ /^(\d+)M$/) {
        $hash{$F[0]}++;
    }
    if ($line =~ /\s+NM:i:(\d+)\s+/) {
        $distance = $1;
        if ($distance > $max_mismatch) {
            $hash{$F[0]}++;
        }
    }
}
close $bam_in;

open $bam_in, "-|", "$samtools view $bam_file" or die "Could not open BAM file: $!\n";
open my $sam_out, ">", "$sample_prefix.filtered.sam" or die "Could not open output SAM file: $!\n";
while (<$bam_in>) {
    chomp;
    my $line = $_;
    my @F = split /\s+/, $line;

    # Apply filtering conditions
    next if ($F[6] ne "=");                  # Only keep pairs mapped to the same chromosome
    next if ($F[4] < $min_map_qual);         # Mapping quality cutoff
    next if ($F[1] & 0x0400);                # Exclude PCR/optical duplicates
    next if ($F[1] & 0x4);                   # Exclude unmapped segments
    next if ($F[1] & 0x8);                   # Exclude pairs with unmapped segments
    next if ($F[1] & 0x10 && $F[1] & 0x20);  # Exclude if both reads are reverse complemented
    next unless ($F[1] & 0x10 || $F[1] & 0x20); # Keep if only one read is reverse complemented
    next if ($F[1] & 0x100);                 # Exclude secondary alignments
    next if ($F[8] > 600 || $F[8] < -600);   # Exclude fragments >600 bp in length
    next if (exists $hash{$F[0]} || exists $hash_pair{$F[0]});          # Exclude fragments with multiple mismatches

    # Final check and output clean SAM lines
    print $sam_out "$line\n";
}

close $bam_in;
close $sam_out;

# Step 2: Calculate fragment lengths for filtered reads from SAM file and output in BED format
open my $bed_out, ">", "$sample_prefix.filtered.bed" or die "Could not open output bed file: $!\n";
open my $sam_in, "<", "$sample_prefix.filtered.sam" or die "Could not open filtered SAM file: $!\n";

while (<$sam_in>) {
    chomp;
    my @F = split /\s+/;
    next if ($F[6] ne "=");
    if ($F[8] < 0 && $F[3] >= $F[7]) {
        my $s = min($F[3], $F[7]) - 1;
        my $len = abs($F[3] - $F[7]) + length($F[9]);
        my $e = $s + $len;
        print $bed_out join("\t", $F[2], $s, $e, $len), "\n";
    }
}

close $sam_in;
close $bed_out;

# Step 3: Extract motif sequences from reference genome

my %hg;
open my $HG, '<', $fa or die "Could not open genome file '$fa': $!\n";
my $S = "";

while (<$HG>) {
    chomp;
    if (/^>/) {
        $_ =~ s/>//g;
        my @fields = split;
        $S = $fields[0];
        $hg{$S} = "";
    } else {
        $hg{$S} .= uc $_;
    }
}
close $HG;

open my $bed_in, '<', "$sample_prefix.filtered.bed" or die "Could not open filtered bed file: $!\n";
open my $motif_out, ">", "$sample_prefix.motif.bed" or die "Could not open motif output file: $!\n";

my $for_count=0;
my $rev_count=0;

while (<$bed_in>) {
    chomp;
    my $line=$_;
    my @F=split/\s+/,$line;
    my $start_forward = $F[1] - $motif_length;
    next if ($F[0] eq "chrM" && $start_forward < 0);
    my $E5 = substr($hg{$F[0]}, $F[1], $motif_length);
    my $x = $F[2] - $motif_length;
    my $E3 = substr($hg{$F[0]}, $x, $motif_length);
    my $com_seq = $E5 . $E3;
    if ($com_seq!~/N/){
      print $motif_out join("\t", $line, $E5, $E3)."\n";
      }
}
close $bed_in;
close $motif_out;

# Step 4: Calculate motif frequencies
my %hash_motif;
my @MER;
open my $motif_in, '<', $motif_list or die "Could not open motif list: $!\n";
while (<$motif_in>) {
    chomp;
    next if (/Type/);
    push @MER, $_;
}
close $motif_in;

my $sum_site = 0;
open my $motif_file, '<', "$sample_prefix.motif.bed" or die "Could not open motif bed file: $!\n";
while (<$motif_file>) {
    chomp;
    my @F = split;
    $sum_site += 2;
    my $motif_up = substr($F[4], 0, 4);
    my $motif_down=revcom($F[5]);
    $motif_down=substr($motif_down,0,4);
    $hash_motif{$motif_up}++;
    $hash_motif{$motif_down}++;

}

close $motif_file;

open my $motif_freq_out, '>', "$sample_prefix.4mer.motif.freq" or die "Could not open motif frequency file: $!\n";
print $motif_freq_out "Type\t$sample_prefix\n";
foreach my $i(@MER){
    if (exists $hash_motif{$i}){
            my $ratio=$hash_motif{$i}/$sum_site;
            print $motif_freq_out "$i\t$ratio\n";
                                        }
        else{
            print $motif_freq_out "$i\t0\n";
        }
}
close $motif_freq_out;

# Step 5: Calculate motif diversity score
open my $freq_in, '<', "$sample_prefix.4mer.motif.freq" or die "Could not open motif frequency file: $!\n";
open my $mds_out, '>', "$sample_prefix.4mer.motif.MDS" or die "Could not open motif frequency file for diversity: $!\n";

my $divScore = 0;
while (<$freq_in>) {
    chomp;
    next if (/Type/);
    my $line=$_;
   my @F=split/\s+/;
   if ($F[1]==0){
	$F[1]=10**-10;
   }
   my $motif_score=-$F[1]*log($F[1])/log(256);
   $divScore+=$motif_score;
}
close $freq_in;

print $mds_out "Sample\t$sample_prefix\tDiversity Score: $divScore\n";


# Helper function for reverse complement
sub revcom {
    my ($seq) = @_;
    $seq =~ tr/ACGTacgt/TGCAtgca/;
    return scalar reverse($seq);
}

```

