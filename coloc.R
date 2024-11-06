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