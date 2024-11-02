library(dplyr)
library(stringr)
library(MatrixEQTL)
library(data.table)
library(argparse)

##########################################################################

parser <- ArgumentParser()

parser$add_argument('--snploc')
parser$add_argument('--snpMat2')
parser$add_argument('--snpMatT')
parser$add_argument('--expMat')
parser$add_argument('--exp_path')
parser$add_argument('--out')
parser$add_argument('--genes_bed')
parser$add_argument('--covars')
trans<- parser$add_mutually_exclusive_group()
trans$add_argument('--trans', action='store_true')
trans$add_argument('--no_trans',action='store_false')

args <- parser$parse_args()

########################################################################

snpMatT<- read.delim(args$snpMatT,header=TRUE,sep=',')
snpMatrix2<- read.delim(args$snpMat2,header=TRUE,sep=',')
exp_tis_3 <- read.table(args$expMat,header=TRUE,stringsAsFactors=FALSE)
tis_expression <- read.delim(args$exp_path, header=TRUE, comment.char="",skip = 2 )   

dir.create(args$out,mode="0777")

write.table(snpMatrix2,paste0(args$out,"/snpMatrix.txt"),row.names = FALSE,quote = FALSE) 

############# Gene location ############################################
Genes_bed <- read.csv(args$genes_bed)
tis_bed<-select(tis_expression,c(Name,Description))
tis_bed<-merge(x=tis_bed,y=Genes_bed,by.x="Name",by.y="gene_id")
tis_bed<-select(tis_bed,!c(Name,X))

tis_bed2<-tis_bed[match(exp_tis_3$Gene,tis_bed$Description),]

########### Create snpMatT #############################################

exp_tis_T<-t(exp_tis_3)
colnames(exp_tis_T)<-exp_tis_3$Gene
exp_tis_T<-exp_tis_T[-1,] 
exp_tis_T<-data.frame(exp_tis_T)

snpMatT2<-snpMatT[,-1]
rownames(snpMatT2)<-snpMatT[,1]


tMatrix_chr<-merge(snpMatT2,exp_tis_T,by="row.names")



write.table(tMatrix_chr,paste0(args$out,"/tMatrix.txt"),quote = FALSE,row.names = FALSE,sep = ",")


################################################
################################################
####### Start MATRIXEQTL pipeline ##############
################################################

SNP_file_name = paste0(args$out,"/snpMatrix.txt")


### Load files

snps = SlicedData$new()
snps$fileDelimiter = " "
snps$fileOmitCharacters= "NA"
snps$fileSkipRows=1
snps$fileSkipColumns=1
snps$fileSliceSize= 2000
snps$LoadFile(SNP_file_name)



gene = SlicedData$new()
gene$fileDelimiter = " "
#gene$fileDelimiter = "\t"
gene$fileOmitCharacters= "NA"
gene$fileSkipRows=1
gene$fileSkipColumns=1
gene$fileSliceSize= 2000
gene$LoadFile(args$expMat)



cvrt = SlicedData$new()
cvrt$fileDelimiter = " "
#cvrt$fileDelimiter = "\t"
cvrt$fileOmitCharacters= "NA"
cvrt$fileSkipRows=1
cvrt$fileSkipColumns=1
if (length(args$covars)>0){
  cvrt$LoadFile(args$covars)
}

###### Check for genes in chromosome for cis ############

  ### Select model type
  
useModel = modelLINEAR
  
  
output_file_name_cis = tempfile()
output_file_name_tra = tempfile()
  
  
  ##### Set p-value threshold, depending on the dataset
pvOutputThreshold_cis = 0.1


if (args$trans){
  pvOutputThreshold_tra = 1e-3
} else {
  pvOutputThreshold_tra = 0 
}



  
errorCovariance = numeric()
  
cisDist= 1e6
  
  
snploc= read.table(args$snploc,header = TRUE, stringsAsFactors = FALSE)
  
  ############################################################################# 
  #############################################################################
  ######### Evaluate chromosome, if it has cis eGene ##########################
  #############################################################################
  
me = Matrix_eQTL_main(
  snps = snps,
  gene = gene,
  cvrt = cvrt,
  output_file_name     = output_file_name_tra,
  pvOutputThreshold     = pvOutputThreshold_tra,
  useModel = useModel,
  errorCovariance = errorCovariance,
  verbose = TRUE,
  output_file_name.cis = output_file_name_cis,
  pvOutputThreshold.cis = pvOutputThreshold_cis,
  snpspos = snploc,
  genepos = tis_bed2,
  cisDist = cisDist,
  pvalue.hist = TRUE,
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE)
  
unlink(output_file_name_tra)
unlink(output_file_name_cis)
  
  
cisQT_chr <- data.frame(me$cis$eqtls)

if (args$trans){
  transQT_chr <- data.frame(me$trans$eqtls)
} 
  
  ################ add for cis ########################
  
por<- c()

  
  
  
for (var in cisQT_chr$snps) {
  
  p<-table(tMatrix_chr[var])
  
  if(length(p)>3){
    p<-p[-4]
  }
  p<-min(p)
  por<-append(por,p)
  
}
  


cisQT_chr["minv"]<-por
  
####### add for trans #########
if (args$trans){
  
  pors<- c()
  
  for (var in transQT_chr$snps) {
    
    p<-table(tMatrix_chr[var])
    
    if(length(p)>3){
      p<-p[-4]
    }
    p<-min(p)
    pors<-append(pors,p)
    
  }
  
  transQT_chr["minv"]<-pors
} 

  

cisQT_chr$beta_se<- cisQT_chr$beta/cisQT_chr$statistic
write.table(cisQT_chr,paste0(args$out,"/cisQT.txt"),quote = FALSE,row.names = FALSE)


if (args$trans){
  
  transQT_chr$beta_se<- transQT_chr$beta/transQT_chr$statistic
  write.table(transQT_chr,paste0(args$out,"/transQT.txt"),quote = FALSE,row.names = FALSE)

}


