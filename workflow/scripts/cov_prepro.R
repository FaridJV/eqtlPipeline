

#############################################################
#############################################################
############################################################
############################################################

library(argparse) #opt_args
library(dplyr)

parser <- ArgumentParser()
parser$add_argument("--inpath")
parser$add_argument("--outpath")
parser$add_argument("--age")
args <- parser$parse_args()


tis_covariates <- read.delim(args$inpath)

p<-rowSums(tail(tis_covariates,n=3)[,-1])
q<-length(colnames(tis_covariates))-1
if (p[1]==q | p[1]==0 ){
  nam<-as.numeric(names(p[1]))
  tis_covariates<-tis_covariates[-c(nam),]
}

if (p[2]==q |p[2]==0  ){
  nam<-as.numeric(names(p[2]))
  tis_covariates<-tis_covariates[-c(nam),]
}
if (p[3]==q*2 | p[3]==q){
  nam<-as.numeric(names(p[3]))
  tis_covariates<-tis_covariates[-c(nam),]
}


age_data <- read.delim(args$age,row.names = 1)
#age_data <- read.delim("D:/GTEX4/Covariates/age_data.txt",row.names = 1)
age_data <- t(age_data)
age_data<- data.frame(age_data)
colu<-colnames(tis_covariates)[-1]
age_data<- age_data %>% select(all_of(colu))
colu<-c("covars",colu)

colnames(tis_covariates)<-colu
age_data<-cbind(data.frame(covars="age_data"),age_data)


tis_covariates2<-rbind(tis_covariates,age_data)
#print(length(tis_covariates2$covars))
write.table(tis_covariates2,args$outpath,row.names = FALSE,quote = FALSE)








