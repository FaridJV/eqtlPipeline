#### Import libraries ####
#options(error=traceback)
library(dplyr)
library(stringr)
library(MatrixEQTL)
library(data.table)
library(reticulate)
library(argparse)

##############################
##### inputs #########

parser <- ArgumentParser()

parser$add_argument('--exp_path')
parser$add_argument('--vcf_path')
parser$add_argument('--out')
parser$add_argument('--genes_bed')
parser$add_argument('--covars')
trans<- parser$add_mutually_exclusive_group()
trans$add_argument('--trans', action='store_true')
trans$add_argument('--no_trans',action='store_false')

args <- parser$parse_args()



########### Import data sets ######
##################################
### Import expression file #######


tis_expression <- read.delim(args$exp_path, header=TRUE, comment.char="",skip = 2 )      ## specify tissue
                 
# Keep only gene and adjust header#

exp_tis_<-select(tis_expression,!c(id,Name))
head<-colnames(exp_tis_)

alg2<-function(x){
  if(str_detect(x,"chr")[1]){ 
    return(x)
  }else{
    
    return(paste("GTEX",unlist(strsplit(x,split = '.',fixed = TRUE))[2],sep = "."))
    
  }}

head<-lapply(head[-1],alg2)
head<-append(head,"Gene",0)

colnames(exp_tis_)<-head

###################

exp_tis_ <- exp_tis_[rowSums(exp_tis_[,-1])>(0.1*(ncol(exp_tis_)-1)),] # Filter genes that have no expression
exp_tis_3<-log(exp_tis_[,-1] + 1) # log transformation
Gene<-exp_tis_$Gene
exp_tis_3<-cbind(Gene,exp_tis_3)

#####################

##### Import previous gene expression files #


##### Import genotype file #######


chr_vcf<-read.delim(args$vcf_path,header = TRUE,comment.char = "") 

snpMatrix<-select(chr_vcf,!c(X.CHROM,POS,REF,ALT,QUAL,FILTER,INFO,FORMAT))
snploc<-select(chr_vcf,c(ID,X.CHROM,POS))

alg<-function(x){
  if(str_detect(x,"chr")[1]){ 
    return(x)
  }else{
    
    return((substr(x,1,3)))
  }}

snpMatrix<-sapply(snpMatrix,alg)

gecode<-c("0/0"="0","1/0"="1","0/1"="1","1/1"="2","./."="NA")

snpMatrix <- data.frame(snpMatrix)
snpMatrix <- snpMatrix %>% mutate_all(funs(str_replace_all(.,gecode)))
samples2<-colnames(snpMatrix)
samples2<-samples2[-1]

##### Check for samples also in expression data #####

samples<-colnames(exp_tis_) 
samples<-samples[-1]
finalsamples<-intersect(samples,samples2)
final<-append(finalsamples,"ID",0)

snpMatrix2<-snpMatrix %>% select(all_of(final))

###################################

final2<-append(finalsamples,"Gene",0)

exp_tis_3<-exp_tis_3%>% select(all_of(final2))
######################################

#### export files snpMatrix and snplocation #####

###################################################

dir.create(args$out,mode="0777")

write.table(exp_tis_3,paste0(args$out,"/exp.txt"),row.names = FALSE,quote = FALSE)

Genes_bed <- read.csv(args$genes_bed)
tis_bed<-select(tis_expression,c(Name,Description))
tis_bed<-merge(x=tis_bed,y=Genes_bed,by.x="Name",by.y="gene_id")
tis_bed<-select(tis_bed,!c(Name,X))

tis_bed2<-tis_bed[match(exp_tis_3$Gene,tis_bed$Description),]



remove(Genes_bed)

################################################



remove(chr_vcf)
remove(exp_tis_)
remove(snpMatrix)
####################################################
####################################################
####################################################
######                             #################
######   Filter variants in python #################
######                             #################
####################################################
####################################################
####################################################


pd<-import("pandas")
snpMat=r_to_py(snpMatrix2)


py_run_string("snpMat=r.snpMat.set_index('ID')")

remove(snpMat)

py_run_string("snpMaT=snpMat.transpose()")

py_run_string("var1=[]")
py_run_string("snpMaT2=[]")

py_run_string(
  "for i in range(snpMaT.shape[1]):
    if min(snpMaT.iloc[:,i].value_counts())<3 or not(set(['0','1','2']).issubset(snpMaT.iloc[:,i].value_counts().axes[0].tolist())) :
      var1.append(i)
")
py_run_string("snpMaT2=snpMaT.drop(snpMaT.columns[var1],axis=1)")
py_run_string("snpMat=snpMaT2.T")


snpMatT<-py$snpMaT2
snpMatrix2<-py$snpMat

py_run_string("del snpMat")
py_run_string("del snpMaT2")
py_run_string("del snpMaT")

snpMatrix2<-cbind(rownames(snpMatrix2),data.frame(snpMatrix2,row.names=NULL))



write.table(snpMatrix2,paste0(args$out,"/snpMatrix.txt"),row.names = FALSE,quote = FALSE)          #specify cromosome


########### Import snpMatT #########################

exp_tis_T<-t(exp_tis_3)
colnames(exp_tis_T)<-exp_tis_3$Gene
exp_tis_T<-exp_tis_T[-1,] #exp_tis_T<- exp_tis_T %>% slice(-1)
exp_tis_T<-data.frame(exp_tis_T)


tMatrix_chr<-merge(snpMatT,exp_tis_T,by="row.names")



write.table(tMatrix_chr,paste0(args$out,"/tMatrix.txt"),quote = FALSE,row.names = FALSE,sep = ",")



################################################
################################################
####### Start MATRIXEQTL pipeline ##############
################################################

SNP_file_name = paste0(args$out,"/snpMatrix.txt")
#paste(paste0("D:/GTEX4/snpMatrix/snpMatrix_chr",chr,"_",tis, ".txt"),sep=" ")
expression_file_name = paste0(args$out,"/exp.txt")
#paste(paste0("D:/GTEX4/Expression2/exp_",tis, "_3.txt"),sep = " ")


#### If no covariates set to character ()


#covararg <- grep("--covars", ca)
#covariates_file_name = ifelse(length(covararg) > 0, ca[covararg+1], "")

  


#paste(paste0("D:/GTEX4/Covariates/",tis, "_covariates.txt"),sep = " ")

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
gene$LoadFile(expression_file_name)



cvrt = SlicedData$new()
#cvrt$fileDelimiter = " "
cvrt$fileDelimiter = "\t"
cvrt$fileOmitCharacters= "NA"
cvrt$fileSkipRows=1
cvrt$fileSkipColumns=1
if (length(covariates_file_name)>0){
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
  pvOutputThreshold_tra = 1e-6
} else {
  pvOutputThreshold_tra = 0 
}



  
errorCovariance = numeric()
  
cisDist= 1e6
  
  
#snpspos= read.table(snps_location_file_name,header = TRUE, stringsAsFactors = FALSE)
#genepos= read.table(gene_location_file_name,header = TRUE, stringsAsFactors = FALSE)
  
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
  transsQT_chr <- data.frame(me$trans$eqtls)
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



