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
#parser.add_argument('--trans',action=argparse.BooleanOptionalAction)
parser.add_argument('--transoutpath')

args  = parser.parse_args() 

#### Import files ########

tMatrix   = pd.read_csv(args.inpath+'/tMatrix.txt',sep=',',header=0,low_memory=False)

### For cis #####

cisQT_chr = pd.read_csv(args.inpath+'/cisQT.txt',sep=" ")
pasz_cis  = []

genes_cis = cisQT_chr['gene'].str.replace('-','.',regex=False)
snps_cis  = cisQT_chr['snps']



for i in range(len(snps_cis)):
    
  q=tMatrix.groupby([snps_cis[i]])[genes_cis[i]].median()
  q1=np.argsort(list(q))
    
  if (q1[1]!=1):
    pasz_cis.append(i)
    

#cisQT_chr = cisQT_chr.drop(pasz_cis)
cisQT2_chr = cisQT_chr.drop(pasz_cis)
cisQT2_chr["linear"]=["yes"]*len(cisQT2_chr)
cisQT3_chr = cisQT_chr[cisQT_chr.index.isin(pasz_cis)]
cisQT3_chr["linear"]=["no"]*len(cisQT3_chr)
cisQT_chr = pd.concat([cisQT2_chr,cisQT3_chr],ignore_index = True)
### output ########
cisQT_chr.to_csv(args.outpath,index=False)


### Check for trans  ##########

if (args.transoutpath):
    
  ### For Trans- #########

  transQT_chr = pd.read_csv(args.inpath+'/transQT.txt',sep=" ")
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

  #transQT_chr = transQT_chr.drop(pasz_trans)
  transQT2_chr = transQT_chr.drop(pasz_trans)
  transQT2_chr["linear"]=["yes"]*len(transQT2_chr)
  transQT3_chr = transQT_chr[transQT_chr.index.isin(pasz_trans)]
  transQT3_chr["linear"]=["no"]*len(transQT3_chr)
  transQT_chr = pd.concat([transQT2_chr,transQT3_chr],ignore_index = True)
  transQT_chr.to_csv(args.transoutpath,index=False)
  



  




