


import pandas as pd
import argparse

parser = argparse.ArgumentParser()

parser.add_argument('--inpath')
parser.add_argument('--outsnpmaT2')
parser.add_argument('--outsnpmat')

args  = parser.parse_args() 

#snpMat=r_to_py(snpMatrix2)
snpMat = pd.read_csv(args.inpath,sep=" ",dtype=str)
snpMat = snpMat.set_index('ID')
#py_run_string("snpMat=r.snpMat.set_index('ID')")

snpMaT=snpMat.transpose()

var1=[]
snpMaT2=[]


for i in range(snpMaT.shape[1]):
    if min(snpMaT.iloc[:,i].value_counts())<3 or not(set(['0','1','2']).issubset(snpMaT.iloc[:,i].value_counts().axes[0].tolist())) :
        var1.append(i)

snpMaT2=snpMaT.drop(snpMaT.columns[var1],axis=1)
snpMat=snpMaT2.T


#snpMatT<-py$snpMaT2
#snpMatrix2<-py$snpMat

snpMaT2.to_csv(args.outsnpmaT2)
snpMat.to_csv(args.outsnpmat)

#snpMatrix2<-cbind(rownames(snpMatrix2),data.frame(snpMatrix2,row.names=NULL))
#write.table(snpMatrix2,paste0(args$out,"/snpMatrix.txt"),row.names = FALSE,quote = FALSE)          #specify cromosome
