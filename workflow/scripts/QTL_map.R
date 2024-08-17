#### Import libraries ####

library(tidyverse)
#library(ggplot2)
library(stringr)
#library(gridExtra)
library(MatrixEQTL)
library(data.table)
library(reticulate)


############################
ca <- commandArgs()

##### inputs #########

#tisarg <- grep("--tissue", ca)
#tis <- ifelse(length(tisarg) > 0, ca[tisarg+1], "Pituitary")

#chrarg <- grep("--chrom", ca)
#chr <- ifelse(length(chrarg) > 0, ca[chrarg+1], "22")

tis_patharg <- grep("--exp-path", ca)
tis_path<- ifelse(length(tis_patharg) > 0, ca[tis_patharg+1], "")

vcf_patharg <- grep("--vcf-path", ca)
vcfpath<-ifelse(length(vcf_patharg) > 0, ca[vcf_patharg+1], "")

pathtemparg <- grep("--out", ca)
pathtemp<-ifelse(length(pathtemparg) > 0, ca[pathtemparg+1], "")

genesbedarg <- grep("--genes-bed", ca)
Genes_bed<-ifelse(length(genesbedarg) > 0, ca[genesbedarg+1], "")
########### Import data sets ######
##################################
### Import expression file #######


tis_expression <- read.delim(tis_path, header=TRUE, comment.char="",skip = 2 )      ## specify tissue

#tis_expression <- read.delim(paste0("D:/GTEX4/Expression2/gene_tpm_2017-06-05_v8_",tis, ".gct.gz"), header=TRUE, comment.char="",skip = 2 )      ## specify tissue
                  
# Keep only gene and adjust header#

exp_tis_<-select(tis_expression,!c(id,Name))
head<-colnames(exp_tis_)

alg2<-function(x){
  if(str_detect(x,"chr")[1]){ 
    return(x)
  }else{
    
    return(paste("GTEX",unlist(strsplit(x,split = '.',fixed = TRUE))[2],sep = "."))
    #return((substr(x,1,10)))
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


chr_vcf<-read.delim(vcfpath,header = TRUE,comment.char = "",skip = 82) 

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

#write.table(snploc,paste0("D:/GTEX4/snpMatrix/snploc_chr",chr,"_",tis, ".txt"),row.names = FALSE,quote = FALSE)

###################################################

write.table(exp_tis_3,paste0(pathtemp,"/exp.txt"),row.names = FALSE,quote = FALSE)

Genes_bed <- read.csv(Genes_bed)
tis_bed<-select(tis_expression,c(Name,Description))
tis_bed<-merge(x=tis_bed,y=Genes_bed,by.x="Name",by.y="gene_id")
tis_bed<-select(tis_bed,!c(Name,X))
#tis_bed2<-merge(x=tis_bed,y=exp_tis_3$Gene,by.x="Description",by.y="Gene")
tis_bed2<-tis_bed[match(exp_tis_3$Gene,tis_bed$Description),]


#write.table(tis_bed2,paste0("D:/GTEX4/Expression2/",tis,"_bed2.txt"),row.names = FALSE,quote = FALSE)

#remove(tis_bed2)
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

#py_run_string("snpMat=r.snpMatrix2.set_index('ID')")
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



write.table(snpMatrix2,paste0(pathtemp,"/snpMatrix.txt"),row.names = FALSE,quote = FALSE)          #specify cromosome


########### Import snpMatT #########################

exp_tis_T<-t(exp_tis_3)
colnames(exp_tis_T)<-exp_tis_3$Gene
exp_tis_T<-exp_tis_T[-1,] #exp_tis_T<- exp_tis_T %>% slice(-1)
exp_tis_T<-data.frame(exp_tis_T)


tMatrix_chr<-merge(snpMatT,exp_tis_T,by="row.names")



write.table(tMatrix_chr,paste0(pathtemp,"/tMatrix.txt"),quote = FALSE,row.names = FALSE,sep = ",")



################################################
################################################
####### Start MATRIXEQTL pipeline ##############
################################################

SNP_file_name = paste0(pathtemp,"/snpMatrix.txt")
#paste(paste0("D:/GTEX4/snpMatrix/snpMatrix_chr",chr,"_",tis, ".txt"),sep=" ")
expression_file_name = paste0(pathtemp,"/exp.txt")
#paste(paste0("D:/GTEX4/Expression2/exp_",tis, "_3.txt"),sep = " ")


#### If no covariates set to character ()
covararg <- grep("--covars", ca)
covariates_file_name = ifelse(length(covararg) > 0, ca[covararg+1], "")

  


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
  cvrt$LoadFile(covariates_file_name)
}

###### Check for genes in chromosome for cis ############

#if (chr %in% chr_cis){
  
#gene_location_file_name = paste(paste0("D:/GTEX4/Expression2/",tis, "_bed2.txt"),sep = "")
 
#snps_location_file_name = paste(paste0("D:/GTEX4/snpMatrix/snploc_chr",chr,"_",tis, ".txt"),sep="") 
  
  ### Select model type
  
useModel = modelLINEAR
  
  
output_file_name_cis = tempfile()
output_file_name_tra = tempfile()
  
  
  ##### Set p-value threshold, depending on the dataset
pvOutputThreshold_cis = 1e-3
pvOutputThreshold_tra = 0 #1e-5
  
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
  
#plot(me)
  
cisQT_chr <- data.frame(me$cis$eqtls)
  #transQT_chr <-data.frame(me$trans$eqtls)
  
py_run_string("import pandas as pd")
py_run_string("import numpy as np")
  
  ################# Import data tMatrix with expression ##################

py_run_string("tMatrix=pd.read_csv(r.pathtemp + '/tMatrix.txt',sep=',',header=0,low_memory=False)")

#py_run_string("tMatrix=pd.read_csv('D:/GTEX4/Expression2/tMatrix_chr'+r.chr + '_' + r.tis + '.txt',sep=',',header=0,low_memory=False)")
  
  ################### Filter for cis #############  
  
py_run_string("pasz_cis=[]")
  
py_run_string("genes_cis=r.cisQT_chr['gene']")
py_run_string("snps_cis=r.cisQT_chr['snps']")
  
py_run_string("for i in range(len(genes_cis)):
  if '-' in genes_cis[i]:
              genes_cis[i]=genes_cis[i].replace('-','.')")

py_run_string("for i in range(len(snps_cis)):
  q=tMatrix.groupby([snps_cis[i]])[genes_cis[i]].median()
  q1=np.argsort(list(q))
  if (q1[1]!=1):
    pasz_cis.append(i)")
py_run_string("r.cisQT_chr = r.cisQT_chr.drop(pasz_cis)")
  
  ################### Filter for trans #############
  #py_run_string("pasz_trans=[]")
  #py_run_string("genes_trans=r.transQT_chr['gene']")
  #py_run_string("snps_trans=r.transQT_chr['snps']")
  #py_run_string("for i in range(len(snps_trans)):
  #q=tMatrix.groupby([snps_trans[i]])[genes_trans[i]].median()
  #q1=np.argsort(list(q))
  #if (q1[1]!=1):
  #  pasz_trans.append(i)")
  #py_run_string("r.transQT_chr = r.transQT_chr.drop(pasz_trans)")
  
  #py_run_string("del tMatrix")
  
  #################################################################
  ##### Add min samples column and export txt file ################
  #################################################################
  
  ################# add for trans ##################
  #por<- c()
  
  #for (var in transQT_chr$snps) {
  #  p<-table(tMatrix_chr[var])
  #  if(length(p)>3){
  #    p<-p[-4]
  #  }
  #  p<-min(p)
  #  por<-append(por,p)
  #}

  #transQT_chr["minv"]<-por
  
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
  
  
  #write.table(transQT_chr,paste0("D:/GTEX4/eqtls/transQT",chr,"_",tis, ".txt"),quote = FALSE,row.names = FALSE)

cisQT_chr$beta_se<- cisQT_chr$beta/cisQT_chr$statistic
write.table(cisQT_chr,paste0(pathtemp,"/cisQT.txt"),quote = FALSE,row.names = FALSE)
  
  
  #} else {

############################################################################# 
#############################################################################
######### Evaluate chromosome, if it has only trans eGenes ##################
#############################################################################

#output_file_name = tempfile()
#me2=Matrix_eQTL_engine(snps, 
#                   gene, 
#                   cvrt = cvrt, 
#                   output_file_name, 
#                   pvOutputThreshold = 1e-5, 
#                   useModel = modelLINEAR, 
#                   errorCovariance = numeric(), 
#                   verbose = TRUE,
#                   pvalue.hist = TRUE,
#                   min.pv.by.genesnp = FALSE,
#                   noFDRsaveMemory = FALSE)



#unlink(output_file_name)

#plot(me2)
#transQT_chr <-data.frame(me2$all$eqtls)

#### Filter by Medians###########

#py_run_string("import pandas as pd")
#py_run_string("import numpy as np")
#py_run_string("pasz=[]")
#py_run_string("genes=r.transQT_chr['gene']")
#py_run_string("snps=r.transQT_chr['snps']")
#py_run_string("tMatrix=pd.read_csv('D:/GTEX4/Expression2/tMatrix_chr'+r.chr + '_' + r.tis + '.txt',sep=' ',header=0,low_memory=False)")
#py_run_string("for i in range(len(snps)):
#  q=tMatrix.groupby([snps[i]])[genes[i]].median()
#  q1=np.argsort(list(q))
#  if (q1[1]!=1):
#    pasz.append(i)")
#py_run_string("r.transQT_chr = r.transQT_chr.drop(pasz)")
#py_run_string("del tMatrix")


#################################################################
##### Add min samples column and export txt file ################
#################################################################

################# add for trans ##################
#por<- c()

#for (var in transQT_chr$snps) {
#  p<-table(tMatrix_chr[var])
#  if(length(p)>3){
#    p<-p[-4]
#  }
#  p<-min(p)
#  por<-append(por,p)
#}

#transQT_chr["minv"]<-por



#write.table(transQT_chr,paste0("D:/GTEX4/eqtls/transQT",chr,"_",tis, ".txt"),quote = FALSE,row.names = FALSE)

 # }

############# Plot some variants ###################
#ggplot(tMatrix_chr,aes(x= chr19_1185473_G_A_b38,y=ADCY3))+geom_boxplot(aes(colour=chr19_1185473_G_A_b38),show.legend = FALSE)+geom_jitter(shape=1,position = position_jitter(0.2),aes(colour=chr19_1185473_G_A_b38),show.legend = FALSE)
#tMatrix_chr$CLYBL<- as.numeric(tMatrix_chr$CLYBL)
#############################################################################
#############################################################################
############ Save workspace session #########################################
#############################################################################

remove(exp_tis_T)
remove(head)
remove(exp_tis_3)
remove(snploc)
remove(snpMatT)
remove(snpMatrix2)
remove(tis_bed)


#save.image(paste0("D:/GTEX4/chr",chr,"_",tis, "_session.RData"))                                             







