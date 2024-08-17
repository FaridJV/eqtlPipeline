########################
########################
########################
## FILTER BY MEDIANS ###
########################
########################

import pandas as pd
import numpy as np
import argparse

##### Input ######


parser = argparse.ArgumentParser()

parser.add_argument('--inpath')
parser.add_argument('--outpath')
parser.add_argument('--trans',action=argparse.BooleanOptionalAction)

args  = parser.parse_args() 

#chro      = "21"
##tis       = "kidney_cortex"
#npath    = "C:/Users/farid/Documents/R/QTL_scripts/example/"
#outpath   = ""
#trans     = ""

#### Import files ########

#tMatrix   = pd.read_csv(inpath+'tMatrix_chr'+chro + '_' + tis + '.txt',sep=',',header=0,low_memory=False)
tMatrix   = pd.read_csv(args.inpath+'tMatrix_chr.txt',sep=',',header=0,low_memory=False)

### For cis #####

#cisQT_chr = pd.read_csv(inpath+'cisQT'+chro +'_'+ tis +'.txt',sep=" ")
cisQT_chr = pd.read_csv(args.inpath+'cisQT.txt',sep=" ")
pasz_cis  = []

genes_cis = cisQT_chr['gene']
snps_cis  = cisQT_chr['snps']


for i in range(len(genes_cis)):
    
  if '-' in genes_cis[i]:
              genes_cis[i]=genes_cis[i].replace('-','.')
        
for i in range(len(snps_cis)):
    
  q=tMatrix.groupby([snps_cis[i]])[genes_cis[i]].median()
  q1=np.argsort(list(q))
    
  if (q1[1]!=1):
    pasz_cis.append(i)
    

cisQT_chr = cisQT_chr.drop(pasz_cis)

### output ########
cisQT_chr.to_csv(args.outpath + "cisQT.txt")


### Check for trans  ##########

if (args.trans):
    
  ### For Trans- #########

  #transQT_chr = pd.read_csv(inpath+'transQT'+chro +'_'+ tis +'.txt',sep=" ")
  transQT_chr = pd.read_csv(args.inpath+'transQT.txt',sep=" ")
  pasz_trans  = []

  genes_trans = transQT_chr['gene']
  snps_trans  = transQT_chr['snps']

  
  for i in range(len(genes_trans)):
      
    if '-' in genes_trans[i]:
                genes_trans[i]=genes_trans[i].replace('-','.')


  for i in range(len(snps_trans)):
    q=tMatrix.groupby([snps_trans[i]])[genes_trans[i]].median()
    q1=np.argsort(list(q))
    if (q1[1]!=1):
      pasz_trans.append(i)

  transQT_chr = transQT_chr.drop(pasz_trans)
  transQT_chr.to_csv(args.outpath + "transQT.txt")



  



