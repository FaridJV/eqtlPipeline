

#############################################################
#############################################################
############################################################
############################################################

ca <- commandArgs()
tisarg <- grep("--tissue", ca)
chrarg <- ifelse(length(tisarg) > 0, ca[tisarg+1], "22") # chr 22 is the default

inpath   <- grep("--inpath", ca)
inpatharg <- ifelse(length(inpath) > 0, ca[inpath+1], "22") 

outpath   <- grep("--outpath", ca)
outpatharg <- ifelse(length(outpath) > 0, ca[outpath+1], "22") 

library(argparse) #opt_args
library(dplyr)



#tis<- args$tissue   #"ovary"
tis<-chrarg

tis_covariates <- read.delim(paste0(inpatharg,tis,"_covariates.txt"))
#tis_covariates <- read.delim(paste0("D:/GTEX4/Covariates/",tis,"_covariates.txt"))

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


age_data <- read.delim(paste0(inpatharg,"age_data.txt"),row.names = 1)
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
write.table(tis_covariates2,paste0(outpatharg,tis,"_covariates2.txt"),row.names = FALSE,quote = FALSE)








