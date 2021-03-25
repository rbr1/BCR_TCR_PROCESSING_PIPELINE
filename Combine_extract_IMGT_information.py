#!/usr/bin/python
#import math
import sys
import os
import commands
import numpy as np
import pandas as pd


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

def Extract_sequences_for_IMGT(id,dir,batch_name):
  dir_use = dir[0]
  for batch in batch_name: 
    IMGT_file_compressed =dir_use+"ORIENTATED_SEQUENCES/ANNOTATIONS/IMGT_RAW/"+batch+".txz"
    command = "tar Jxvf "+IMGT_file_compressed+" -C "+dir_use+"ORIENTATED_SEQUENCES/ANNOTATIONS/IMGT_RAW/"+batch+"/"
    os.system("mkdir "+dir_use+"ORIENTATED_SEQUENCES/ANNOTATIONS/IMGT_RAW/"+batch+"/")
    os.system(command)
  ids = {}
  print "\n"
  for s in range(len(id)):
    sample,dir_use = id[s],dir[s]
    fasta_file = dir_use+"ORIENTATED_SEQUENCES/NETWORKS/Fully_reduced_"+sample+".fasta"
    fh= open(fasta_file,"r")
    for header,sequence in fasta_iterator(fh):
      ids[header.split("__")[0]] = sample
    print "\r","read samples:",s, "n sequences:",len(ids),
  fh.close()
  dir_IMGT = dir_use+"ORIENTATED_SEQUENCES/ANNOTATIONS/IMGT_SPLIT/"
  ## Added by LEO as caused an error
  os.system("mkdir "+dir_IMGT)
  ##################################
  files = ['10_V-REGION-mutation-hotspots.txt',
  '1_Summary.txt',
  '2_IMGT-gapped-nt-sequences.txt',
  '3_Nt-sequences.txt',
  '4_IMGT-gapped-AA-sequences.txt',
  '5_AA-sequences.txt',
  '6_Junction.txt',
  '7_V-REGION-mutation-and-AA-change-table.txt',
  '8_V-REGION-nt-mutation-statistics.txt',
  '9_V-REGION-AA-change-statistics.txt']
  for f in range(len(files)):
    for s in range(len(id)):
       out_file = dir_IMGT+"IMGT_"+id[s]+"_"+files[f]
       fh=open(out_file, "w")
       fh.close()
    info = {}
    for batch in batch_name:
      input_file = dir_use+"ORIENTATED_SEQUENCES/ANNOTATIONS/IMGT_RAW/"+batch+"/"+files[f]
      fh=open(input_file,"r")
      for l in fh:
        if(l[0]!="S"):
          l1 = l.strip()
          l=l1.split()
          if(l[1].split("__")[0] in ids): 
            sam = ids[l[1].split("__")[0]]
            if(sam in info): info[sam] = info[sam]+l1+"\n"
            else:info[sam]=l1+"\n"
            if(len(info[sam])>10000):
              Write_out(info[sam], dir_IMGT+"IMGT_"+sam+"_"+files[f])
              info[sam] = ''
          else:
            print l
      fh.close()
      for sam in info: 
        Write_out(info[sam], dir_IMGT+"IMGT_"+sam+"_"+files[f])
        info[sam] = ''
    print f, files[f]
  return()

def Write_out(out, file):
  fh=open(file,"a")
  fh.write(out)
  fh.close()
  return()

################################################
args=sys.argv
batch_file = args[1]
id, dir = Get_sample_files(batch_file) # samples.post file as argument 

# Edited by LEO to run on command line without need to edit. 
data = pd.read_csv(batch_file, header = None, delimiter = "\t")
output_dir = data.loc[1, 7]
path_to_IMGTfiles = output_dir + '/ORIENTATED_SEQUENCES/ANNOTATIONS/IMGT_RAW'
batch_name = os.listdir(path_to_IMGTfiles) 
for i in range(len(batch_name)): 
    batch_name[i] = batch_name[i].split('.txz')[0]  

##### extract sequences from batched files for IMGT
Extract_sequences_for_IMGT(id,dir,batch_name)
### check sequences which ones are missed?? 


