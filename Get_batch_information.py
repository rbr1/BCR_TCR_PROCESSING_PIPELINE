#!/usr/bin/python
#import math
import sys
import os
import commands
import numpy as np

def fasta_iterator(fh):
  while True:
    line = fh.readline()
    if line.startswith('>'): break  
  while True:
    header = line[1:-1].rstrip()
    sequence = fh.readline().rstrip()
    while True:
      line = fh.readline()
      if not line: break
      if line.startswith('>'): break
      sequence += line.rstrip()
    yield(header, sequence)
    if not line: return

def Get_sample_files(batch_file):
  id, dir = [],[]
  fh=open(batch_file,"r")
  ind = 0
  for l in fh:
    ind = ind+1
    if(l[0]!="#" and len(l)>3):
      l=l.strip().split()
      id.append(l[0])
      dir.append(l[7])
  fh.close()
  return(id, dir)

def Batch_sequences_for_IMGT(id,dir):
  file_out_pre = dir[0]+"ORIENTATED_SEQUENCES/NETWORKS/Fully_reduced_"+batch_file.replace(".txt","").replace("Samples_","")+"_"
  out, ind = '',0
  batch_size,batch = 1000000-10,1
  n = 0
  fh=open(file_out_pre+str(batch)+".fasta","w")
  fh.close()
  for s in range(len(id)):
    sample,dir_use = id[s],dir[s]
    fasta_file = dir_use+"ORIENTATED_SEQUENCES/NETWORKS/Fully_reduced_"+sample+".fasta"
    fh= open(fasta_file,"r")
    for header,sequence in fasta_iterator(fh):
      n,ind=n+1,ind+1
      out=out+">"+header+"\n"+sequence+"\n"
      if(ind>1000 or n>=batch_size):
        Write_out(out, file_out_pre+str(batch)+".fasta")
        out, ind = '',0
        if(n>=batch_size):
          n,batch = 0,batch+1
          fh=open(file_out_pre+str(batch)+".fasta","w")
          fh.close()
    fh.close()
  Write_out(out, file_out_pre+str(batch)+".fasta")
  out, ind = '',0
  print "scp -p mfj169@rescomp1.well.ox.ac.uk:"+file_out_pre+"* ./"
  return()

def Write_out(out, file):
  fh=open(file,"a")
  fh.write(out)
  fh.close()
  return()

def Get_isotype_depth(id,dir):
  isotypes_uniq,isotypes_total = {},{}
  all_uniq, all_total = {},{}
  for s in range(len(id)):
    sample,dir_use = id[s],dir[s]
    isotype_file = dir_use+"ORIENTATED_SEQUENCES/ANNOTATIONS/IsoTyper_chain_repertoire_statistics_file_"+sample+".txt"
    fh=open(isotype_file,"r")
    for l in fh:
      if(l[0]!="#"):
        l=l.strip().split()
        iso, umis, uniqs = l[1],int(l[2]), int(l[3])
        if(iso in isotypes_uniq):
          isotypes_uniq[iso], isotypes_total[iso] = isotypes_uniq[iso]+[uniqs], isotypes_total[iso] +[umis]
        else:isotypes_uniq[iso], isotypes_total[iso] = [uniqs],[umis]
        if(sample in all_uniq):
          all_uniq[sample],all_total[sample] = all_uniq[sample]+uniqs,all_total[sample]+umis
        else:all_uniq[sample],all_total[sample] = uniqs,umis
    fh.close()
  array_uniq, array_total = [],[]
  for s in all_uniq:
    array_uniq, array_total = array_uniq +[all_uniq[s]], array_total + [all_total[s]]
  isotypes_uniq["all"] = array_uniq
  isotypes_total["all"] = array_total
  out = "#isotype\ttype\tmin\t5th percentile\t10th percentile\t20th percentile\n"
  for iso in isotypes_uniq:
    quantiles = [np.percentile(isotypes_uniq[iso],5), np.percentile(isotypes_uniq[iso],10),np.percentile(isotypes_uniq[iso],20)]
    out = out+"\t".join(map(str, [iso,"UNIQ",min(isotypes_uniq[iso])]+quantiles))+"\n"
  for iso in isotypes_total:
    quantiles = [np.percentile(isotypes_total[iso],5), np.percentile(isotypes_total[iso],10),np.percentile(isotypes_total[iso],20)]
    out = out+"\t".join(map(str, [iso,"TOTAL",min(isotypes_total[iso])]+quantiles))+"\n"
  out_file = dir_use+"ORIENTATED_SEQUENCES/ANNOTATIONS/Sampling_depth_per_isotype_"+batch_file.replace(".txt","")+".txt"
  print out_file
  fh=open(out_file, "w")
  fh.write(out)
  fh.close()
  return()

################################################
args=sys.argv
batch_file = args[1]
print batch_file

id, dir = Get_sample_files(batch_file)


Get_isotype_depth(id,dir)
Batch_sequences_for_IMGT(id,dir)

