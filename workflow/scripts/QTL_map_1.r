#### Import libraries ####
#options(error=traceback)
library(dplyr)
library(stringr)
library(data.table)
library(argparse)

##############################
##### inputs #########

parser <- ArgumentParser()

parser$add_argument('--exp_path')
parser$add_argument('--vcf_path')
parser$add_argument('--outsnpMat')
parser$add_argument('--outexpMat')
parser$add_argument('--outsnploc')


args <- parser$parse_args()



########### Import data sets ######
##################################
### Import expression file #######


tis_expression <- read.delim(args$exp_path, header=TRUE, comment.char="",skip = 2 )      
                 
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

write.table(exp_tis_3,args$outexpMat,row.names = FALSE,quote = FALSE)
write.table(snpMatrix2,args$outsnpMat,row.names = FALSE,quote = FALSE)
write.table(snploc,args$outsnploc,row.names = FALSE,quote = FALSE)