# Generate_repertoire_statistics.py is developed by Rachael Bashford-Rogers (2020)
# at the University of Oxford and Univeristy of Cambridge
# E-mail: rbr1@well.ox.ac.uk

# If you find the methods in ImmuneReceptor_PROCESSING_PIPELINE, please cite the following reference:
# Bashford-Rogers, R. et al. Nature 2019 (https://www.nature.com/articles/s41586-019-1595-3.pdf)

# Copyright (C) 2020  Dr Rachael Bashford-Rogers

#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.

#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.

#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <https://www.gnu.org/licenses/>.

#!/usr/bin/python
import math
import sys
from collections import defaultdict
import os
import re
import time
import commands
import numpy as np
from numpy import outer
from operator import itemgetter, attrgetter, add
import copy
import random
import networkx as nx

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

class Tree(defaultdict):
  def __init__(self, value=None):
    super(Tree, self).__init__(Tree)
    self.value = value

def Uniq(v):
  C=set(v)
  return list(C)

def Get_network_statistics_per_chain(cluster_file, sample, dir,per_chain_repertoire_statistics_file,isotyper_primer_set):
  fh=open(cluster_file,"r")
  cluster,vertices=Tree(), Tree()
  index, totalc, totalv, totalreads, sizesv, c_sizes,vertices_in_max_cluster = 0,0,0,0,[],{},0
  total_v, total_reads = [],[]
  sizesv, c_sizes= {},{}
  chains_short = []
  t1,t2 = 0,0
  n = 0
  for l in fh:
    if(l[0]!="#"):
      l=l.strip().split()
      id = l[2]
      chains, freq, id_short = id.split("|")[1].split("_"), map(int, id.split("__")[1].split("|")[0].split("_")), id.split("__")[0]
      t1 = t1+sum(freq)
      n = n+1
      if(len(chains_short)==0):
        for i in range(0,len(chains)):
          c = chains[i].split("*")[0]
          if(isotyper_primer_set =="INNER_DD"):
            if(c in ["IGHA1","IGHA2"]):c = "IGHA"
            elif(c in ["IGHG1","IGHG2"]):c = "IGHG1/2"
          if(c not in chains_short):chains_short.append(c)
          #if(c in chains_short):chains_index[chains[i].split("*")[0]] = chains_short.index(c)
          #else:
          #  chains_short.append(c)
          #  chains_index[chains[i].split("*")[0]] = chains_short.index(c)
      non_zero = [i for i in range(len(freq)) if freq[i]!=0]
      if(len(total_v)==0):
        total_v, total_reads = [0]*len(chains_short),[0]*len(chains_short)
        for c in chains_short:
          sizesv[c], c_sizes[c] = [],[]
      for i in non_zero:
        c = chains[i].split("*")[0]
        if(isotyper_primer_set =="INNER_DD"):
          if(c in ["IGHA1","IGHA2"]):c = "IGHA"
          elif(c in ["IGHG1","IGHG2"]):c = "IGHG1/2"
        cluster[c][l[1]][freq[i]][id_short].value=1
        index = chains_short.index(c)
        total_v[index] = total_v[index]+1
        sizesv[c] =sizesv[c]+[freq[i]]
        total_reads[index] =total_reads[index]+freq[i]
  fh.close()
  print total_reads, t1, n
  if(t1 != sum(total_reads)):print "ERROR IN COUNTING!!"
  out="#Id\tIsotype\tN reads\tN vertices\tVertex Gini Index\tCluster Gini Index\tLargest Cluster (%)\t2nd Largest Cluster (%)\n"
  #out="#Id\tAnalysis\tN reads\tN vertices\tVertex Gini Index\tCluster Gini Index\tLargest Cluster (%)\t2nd Largest Cluster (%)\t% Vertices in largest cluster\tVertex Renyi\tCluster Renyi\tGene\tSpecies\n"
  for c1 in chains_short: 
    cluster_sizes_sub = []
    for clus in cluster[c1]:
      f = 0
      for f1 in cluster[c1][clus]:
        f = f+(f1*len(cluster[c1][clus][f1]))
      cluster_sizes_sub = cluster_sizes_sub+[f]
    if(len(cluster_sizes_sub)>0):
      (vpoints,vvdf)=VDF(sizesv[c1])
      (cpoints,cvdf)=VDF(cluster_sizes_sub)
      vgini, cgini=Gini_index(cpoints,cvdf, vpoints,vvdf,sum(cluster_sizes_sub),sum(sizesv[c1]),total_v)
      max_pop, max_1_pop = cpoints[len(cpoints)-1]*100.0/sum(sizesv[c1]), cpoints[len(cpoints)-2]*100.0/sum(sizesv[c1])
      out = out+str(sample)+"\t"+c1+"\t"+str(sum(sizesv[c1]))+"\t"+str(len(sizesv[c1]))+"\t"+str(vgini)+"\t"+str(cgini)+"\t"+str(max_pop)+"\t"+str(max_1_pop)+"\n"
  fh=open(per_chain_repertoire_statistics_file, "w")
  fh.write(out)
  fh.close()
  return()

def VDF (n):
  points=sorted(Uniq(n))
  vdf=[]
  for i in range(0,len(points)):
    vdf.append(n.count(points[i]))
  return (points,vdf)

def Get_network_statistics(cluster_file,id,dir,network_statistics,gene,species,cluster_size_distribution,vertex_size_distribution):
  (cpoints,cvdf, vpoints,vvdf,totalc,totalv,totalreads,c_sizes,vertices_in_max_cluster)=Get_cluster_vertex_distributions(cluster_file)
  Print_distributions( vpoints,vvdf, vertex_size_distribution)
  Print_distributions(cpoints,cvdf,cluster_size_distribution)
  (vrenyi,crenyi)=Renyi_entropy(cpoints,cvdf, vpoints,vvdf,totalc,totalv,totalreads)
  (vgini, cgini)=Gini_index(cpoints,cvdf, vpoints,vvdf,totalc,totalv,totalreads)
  (max_pop, max_1_pop)=Proportional_measures(c_sizes, totalreads)
  out="#Id\tAnalysis\tN reads\tN vertices\tVertex Gini Index\tCluster Gini Index\tLargest Cluster (%)\t2nd Largest Cluster (%)\t% Vertices in largest cluster\tVertex Renyi\tCluster Renyi\tGene\tSpecies\n"
  out=out+str(id)+"\tOVERALL\t"+str(totalreads)+"\t"+str(totalv)+"\t"+str(vgini)+"\t"+str(cgini)+"\t"+str(max_pop)+"\t"+str(max_1_pop)+"\t"+str(vertices_in_max_cluster*100.0/totalv)+"\t"+str(vrenyi)+"\t"+str(crenyi)+"\t"+gene+"\t"+species+"\n"
  fh=open(network_statistics, "a")
  fh.write(out)
  fh.close()
  return()

def Get_cluster_vertex_distributions(cluster_file):
  fh=open(cluster_file,"r")
  cluster,vertices=Tree(), Tree()
  index, totalc, totalv, totalreads, sizesv, c_sizes,vertices_in_max_cluster = 0,0,0,0,[],{},0
  for l in fh:
    index=index+1
    if (index>1):
      l=l.strip()
      l=l.split()
      cluster[l[1]][l[2]].value=1
      size = int(l[3])
      sizesv.append(size)
      totalv=totalv+1
      totalreads=totalreads+size
      if(int(l[1])==1):vertices_in_max_cluster=vertices_in_max_cluster+1
      if(l[1] in c_sizes):c_sizes[l[1]]=c_sizes[l[1]]+size
      else:c_sizes[l[1]]=size
  fh.close()
  sizes=[] 
  totalc=len(cluster)
  for c in cluster:
    sizes.append(len(cluster[c]))
  (cpoints,cvdf)=VDF(sizes)
  (vpoints,vvdf)=VDF(sizesv)
  return(cpoints,cvdf, vpoints,vvdf,totalc,totalv,totalreads,c_sizes, vertices_in_max_cluster)
  
def Print_distributions(points,cdf, file_out):
  out = "#size\tfrequency of size\n"
  for i in range(0,len(points)):
    out=out+str(points[i])+"\t"+str(cdf[i])+"\n"
  fh=open(file_out,"w")
  fh.write(out)
  fh.close()
  return()

def Renyi_entropy(cpoints,cvdf, vpoints,vvdf,totalc,totalv,totalreads):
  vrenyi = 0
  crenyi = 0
  tv=totalreads*totalreads*1.0
  tc=totalv*totalv*1.0
  for i in range(0,len(vpoints)):
    vrenyi = vrenyi + vvdf[i]*(vpoints[i]*vpoints[i]/tv)
  for i in range(0,len(cpoints)):
    crenyi = crenyi + cvdf[i]*(cpoints[i]*cpoints[i]/tc)
  return(vrenyi,crenyi)

def Gini_index(cpoints,cvdf, vpoints,vvdf,totalc,totalv,totalreads): 
  (vgini)=Get_Gini(vpoints,vvdf)
  (cgini)=Get_Gini(cpoints,cvdf)
  return(vgini, cgini)

def Proportional_measures(c_sizes, totalreads):
  sizes=[]
  for c in c_sizes:
    sizes.append((c,c_sizes[c]))
  s=sorted(sizes, key=itemgetter(1), reverse=True)
  (max_pop, max_1_pop)=(s[0][1]*100.0/totalreads, s[1][1]*100.0/totalreads)
  return(max_pop, max_1_pop)

def Gini_index(cpoints,cvdf, vpoints,vvdf,totalc,totalv,totalreads): 
  (vgini)=Get_Gini(vpoints,vvdf)
  (cgini)=Get_Gini(cpoints,cvdf)
  return(vgini, cgini)

def Get_Gini(n,v):
  values=[]
  for i in range(0,len(n)):
    for j in range(0,v[i]):
      values.append(n[i])
  n = len(values)
  assert(n > 0), 'Empty list of values'
  sortedValues = sorted(values) #Sort smallest to largest
  cumm = [0]
  for i in range(n):
    cumm.append(sum(sortedValues[0:(i + 1)]))
  LorenzPoints = [[], []]
  sumYs = 0           #Some of all y values
  robinHoodIdx = -1   #Robin Hood index max(x_i, y_i)
  for i in range(1, n + 2):
    x = 100.0 * (i - 1)/n
    y = 100.0 * (cumm[i - 1]/float(cumm[n]))
    LorenzPoints[0].append(x)
    LorenzPoints[1].append(y)
    sumYs += y
    maxX_Y = x - y
    if maxX_Y > robinHoodIdx: robinHoodIdx = maxX_Y   
  giniIdx = 100 + (100 - 2 * sumYs)/n #Gini index 
  return(giniIdx/100)

def Get_constant_region_distribution(seq_file,constant_region_count_file,pat,annot_file):
  fh=open(seq_file,"r")
  #genes = [0]
  #const=[]
  IHGDM_ind,cw_ind = [],[]
  IGHDM_seqs,counts = {},{}
  for header,sequence in fasta_iterator(fh):
    freq, const = map(int,header.split("__")[1].split("|")[0].split("_")), header.split("|")[1].split("_")
    if(len(IHGDM_ind)==0):
      IHGDM_ind = [i for i in range(len(const)) if const[i] in ["IGHD","IGHM"]]
      cw_ind= [i for i in range(len(const)) if const[i] not in ["IGHD","IGHM"]]
      for i in const:
        counts[i] = 0
      counts["class_switched"] = 0
      counts["ALL"] = 0
    nz = [i for i in range(len(freq)) if freq[i]!=0]
    counts["ALL"] = counts["ALL"]+1
    for i in nz: 
      if(i in IHGDM_ind): 
        IGHDM_seqs[header] = 1
      else:counts["class_switched"] = counts["class_switched"]+1
      counts[const[i]] = counts[const[i]]+1
  fh.close()
  counts["IGHD/M_unmutated"] = 0
  counts["IGHD/M_mutated"] = 0
  fh=open(annot_file, "r")
  for l in fh: 
    if(l[0]!="#"):
      l=l.strip().split()
      if(l[0] in IGHDM_seqs):
        if(len(l)>=19):
          mm_v,mm_j =  int(l[17]),int(l[18])
          if(mm_v+mm_j <= 4):mut = "IGHD/M_unmutated"
          else: mut = "IGHD/M_mutated"
          counts[mut] = counts[mut]+1
  fh.close()
  genes,gene_ids = [],[]
  for i in counts: 
    genes,gene_ids = genes+[counts[i]],gene_ids+[i]
  out='#sample\tgene\tfrequency\n'
  if(len(genes)>0):
    for i in range(0,len(genes)):
      out=out+pat+"\t"+gene_ids[i]+"\t"+str(genes[i])+"\n"
  fh=open(constant_region_count_file,"w")
  fh.write(out)
  fh.close()
  return()

def Get_gene_frequencies(annot_file, gene_freq_file,gene,id):
  fh=open(annot_file,"r")
  genes,total,found = {},0,0
  for l in fh:
    if(l[0]!="#"):
      l=l.strip().split()
      if(len(l)>=13):
        if(l[0].count("__")!=0):
          f, v, j = sum(map(int,l[0].split("__")[1].split("|")[0].split("_"))), l[1].split("*")[0], l[13].split("*")[0]
        else:f, v, j = 1, l[1].split("*")[0], l[13].split("*")[0]
        total = total+f
        if(v.count("not_found")==0 and j.count("J")!=0):
          v=v+"|"+j # change back to v+j for complete list
          if(v in genes):genes[v]=genes[v]+f
          else:genes[v]=f
          found = found+f
        else:
          print l
      elif(len(l)==9):
        if(l[0].count("__")!=0):
          f, v, j = sum(map(int,l[0].split("__")[1].split("|")[0].split("_"))), l[1].split("*")[0], l[3].split("*")[0]
        #if(v.count(gene)!=0 and j.count(gene)!=0):
        if(v.count("not_found")==0 and j.count("J")!=0):
          v=v+"|"+j # change back to v+j for complete list
          if(v in genes):genes[v]=genes[v]+f
          else:genes[v]=f
          found = found+f
      else:
        print l
  fh.close()
  print "TOTAL READS:",total,"FOUND READS:", found
  out=''
  for g in genes:
    out=out+id+"\t"+g+"\t"+str(genes[g])+"\t"+g[0:5]+"\n"
  fh=open(gene_freq_file, "w")
  fh.write(out)
  fh.close()
  return()

def Get_cluster_stats(annot_file, cluster_file, id, dir,network_statistics,gene,species,seq_file,loc):
  primer_reference={}
  threshold = 1.0 ### the minimum % of reads in cluster to be reported
  threshold = 0.1
  sig_clust_info,ids,clusters =Get_large_clusters(cluster_file,threshold,100000000000000000000000000000)
  del clusters
  cluster_annotation,CDR3s = Get_annotation_for_clusters(annot_file, ids)
  seqs = Get_max_sequences(seq_file,ids)
  out="#Id\tAnalysis\tCluster ID\t% of Repertoire\tN Reads\tN Vertices\tV gene\tJ gene\tMean mutations in V\tMax Sequence Size (reads)\tMax Sequence ID\tGene\tSpecies\tSpecific Primer\tLargest vertex sequence\n"
  if(len(sig_clust_info)>0):
    for i in range(0,len(sig_clust_info)):
      c,perc_network, number_vertices, number_reads, max_seq, reads_in_max_seq = sig_clust_info[i]
      if(max_seq.split("__")[0] in seqs):seq = seqs[max_seq.split("__")[0]]
      else:seq = max_seq
      if(c in cluster_annotation):max_vj, mean_muts = cluster_annotation[c]
      else:max_vj, mean_muts = "NA\tNA","NA"
      if(max_vj.split("\t")[0] in primer_reference):primer_sequence = primer_reference[max_vj.split("\t")[0]]
      else:primer_sequence="NA"
      out=out+id+"\tCLUSTER_ANALYSIS\t"+str(c)+"\t"+str(perc_network)+"\t"+str(number_reads)+"\t"+str(number_vertices)+"\t"+max_vj+"\t"+str(mean_muts)+"\t"+str(reads_in_max_seq)+"\t"+str(max_seq)+"\t"+gene+"\t"+species+"\t"+primer_sequence+"\t"+seq+"\n"
  else:
    out=out+id+"\tCLUSTER_ANALYSIS\tNO CLUSTERS LARGER THAN "+str(threshold)+"% OF REPERTOIRE\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\t"+gene+"\t"+species+"\tNA\tNA\n"
  fh=open(network_statistics, "a")
  fh.write(out)
  fh.close()
  return()

def Get_max_sequences(seq_file,ids):
  ids_find= {}
  for id in ids:
    ids_find[id.split("__")[0]] = id
  fh=open(seq_file,"r")
  seqs={}
  for header,sequence in fasta_iterator(fh):
    if(header.split("__")[0] in ids_find):
      seqs[header.split("__")[0]]=sequence
  fh.close()
  return(seqs)

def Get_annotation_for_clusters(annot_file, ids):
  ids_find ={}
  for id in ids:
    ids_find[id.split("__")[0]]=id
  fh=open(annot_file,"r")
  cluster_annot,cluster_mutations,CDR3s=Tree(),{},{}
  for l in fh:
    if(l[0]!="#"):
      l=l.strip().split()
      if(l[0].split("__")[0] in ids_find):
        c = ids[ids_find[l[0].split("__")[0]]]
        if(len(l)>=14):
          v,j = l[1],l[13]
          if(j.count("J")==0):print l
          vj= l[1]+"\t"+l[13]
          cluster_annot[c][vj][l[0]].value=1
          if(len(l)>=20 and l[19].count("CDR")==0):
            CDR3s[l[0].split("__")[0]]=l[17]
            if(c in cluster_mutations):cluster_mutations[c] = copy.deepcopy(cluster_mutations[c])+[int(l[19])]
            else:cluster_mutations[c] =[int(l[19])]
  cluster_annotation = {}
  for c in cluster_annot:
    max_vj,max_f_vj = '',0
    for vj in cluster_annot[c]:
      if(len(cluster_annot[c][vj])>max_f_vj):max_vj,max_f_vj = vj, len(cluster_annot[c][vj])
    if(c in cluster_mutations):
      m = cluster_mutations[c]
      mean_muts = sum(m)*1.0/len(m)
    else: mean_muts = "NA"
    cluster_annotation[c] = [max_vj, mean_muts]
  del cluster_annot, cluster_mutations
  return(cluster_annotation,CDR3s)

def Get_large_clusters(cluster_file,threshold,max_number_of_clusters_to_include):
  fh=open(cluster_file,"r")
  total,clusters,clust_size=0,Tree(),{}
  for l in fh:
    if(l[0]!="#"):
      l=l.strip().split()
      c,id,f = l[1], l[2], int(l[3])
      id = id.split("__")[0]+"__"+l[3]
      total = total+f
      if(c in clust_size):clust_size[c]=clust_size[c]+f
      else:clust_size[c]=f
      clusters[c][id][f].value=1
  fh.close()
  sig_clust_info,ids = [],{}
  for c in clust_size:
    prop= clust_size[c]*100.0/total
    if(prop>=threshold):
      perc_network, number_vertices, number_reads, max_seq, reads_in_max_seq = prop, len(clusters[c]),0,'',0
      for id in clusters[c]:
        ids[id]=c
        for f in clusters[c][id]:
          number_reads = number_reads+f
          if(f>reads_in_max_seq):
            max_seq, reads_in_max_seq = id,f
      sig_clust_info.append([c,perc_network, number_vertices, number_reads, max_seq, reads_in_max_seq])
  sig_clust_info =sorted(sig_clust_info,key=itemgetter(1),reverse=True)
  if(len(sig_clust_info)>max_number_of_clusters_to_include):
    ids_new,sig_clust_info_new = {},[]
    for i in range(0,max_number_of_clusters_to_include):
      c = sig_clust_info[i][0]
      for id in clusters[c]:
        ids_new[id]=c
      sig_clust_info_new.append(sig_clust_info[i])
    return(sig_clust_info_new, ids_new, clusters)
  else:
    return(sig_clust_info,ids, clusters)

def Initalise_file(file):
  fh=open(file,"w")
  fh.close()
  return()

def Write_out(out, file):
  fh = open (file,"a")
  fh.write(out)
  fh.close()
  return()

def Intialise_files(dir):
  dirs_to_add=[dir, dir+"TMP/"]
  for d in dirs_to_add:
    c=commands.getoutput("ls "+d)
    if(c.count("No such file or directory")==1):commands.getoutput("mkdir "+d)
  return()

def Get_network_statistics_per_chain(cluster_file, sample, dir,per_chain_repertoire_statistics_file,isotyper_primer_set):
  fh=open(cluster_file,"r")
  cluster,vertices=Tree(), Tree()
  index, totalc, totalv, totalreads, sizesv, c_sizes,vertices_in_max_cluster = 0,0,0,0,[],{},0
  total_v, total_reads = [],[]
  sizesv, c_sizes= {},{}
  chains_short = []
  t1,t2 = 0,0
  n = 0
  for l in fh:
    if(l[0]!="#"):
      l=l.strip().split()
      id = l[2]
      chains, freq, id_short = id.split("|")[1].split("_"), map(int, id.split("__")[1].split("|")[0].split("_")), id.split("__")[0]
      t1 = t1+sum(freq)
      n = n+1
      if(len(chains_short)==0):
        for i in range(0,len(chains)):
          c = chains[i].split("*")[0]
          if(isotyper_primer_set =="INNER_DD"):
            if(c in ["IGHA1","IGHA2"]):c = "IGHA"
            elif(c in ["IGHG1","IGHG2"]):c = "IGHG1/2"
          if(c not in chains_short):chains_short.append(c)
          #if(c in chains_short):chains_index[chains[i].split("*")[0]] = chains_short.index(c)
          #else:
          #  chains_short.append(c)
          #  chains_index[chains[i].split("*")[0]] = chains_short.index(c)
      non_zero = [i for i in range(len(freq)) if freq[i]!=0]
      if(len(total_v)==0):
        total_v, total_reads = [0]*len(chains_short),[0]*len(chains_short)
        for c in chains_short:
          sizesv[c], c_sizes[c] = [],[]
      for i in non_zero:
        c = chains[i].split("*")[0]
        if(isotyper_primer_set =="INNER_DD"):
          if(c in ["IGHA1","IGHA2"]):c = "IGHA"
          elif(c in ["IGHG1","IGHG2"]):c = "IGHG1/2"
        cluster[c][l[1]][freq[i]][id_short].value=1
        index = chains_short.index(c)
        total_v[index] = total_v[index]+1
        sizesv[c] =sizesv[c]+[freq[i]]
        total_reads[index] =total_reads[index]+freq[i]
  fh.close()
  print total_reads, t1, n
  if(t1 != sum(total_reads)):print "ERROR IN COUNTING!!"
  out="#Id\tIsotype\tN reads\tN vertices\tVertex Gini Index\tCluster Gini Index\tLargest Cluster (%)\t2nd Largest Cluster (%)\n"
  #out="#Id\tAnalysis\tN reads\tN vertices\tVertex Gini Index\tCluster Gini Index\tLargest Cluster (%)\t2nd Largest Cluster (%)\t% Vertices in largest cluster\tVertex Renyi\tCluster Renyi\tGene\tSpecies\n"
  for c1 in chains_short: 
    cluster_sizes_sub = []
    for clus in cluster[c1]:
      f = 0
      for f1 in cluster[c1][clus]:
        f = f+(f1*len(cluster[c1][clus][f1]))
      cluster_sizes_sub = cluster_sizes_sub+[f]
    if(len(cluster_sizes_sub)>0):
      (vpoints,vvdf)=VDF(sizesv[c1])
      (cpoints,cvdf)=VDF(cluster_sizes_sub)
      vgini, cgini=Gini_index(cpoints,cvdf, vpoints,vvdf,sum(cluster_sizes_sub),sum(sizesv[c1]),total_v)
      max_pop, max_1_pop = cpoints[len(cpoints)-1]*100.0/sum(sizesv[c1]), cpoints[len(cpoints)-2]*100.0/sum(sizesv[c1])
      out = out+str(sample)+"\t"+c1+"\t"+str(sum(sizesv[c1]))+"\t"+str(len(sizesv[c1]))+"\t"+str(vgini)+"\t"+str(cgini)+"\t"+str(max_pop)+"\t"+str(max_1_pop)+"\n"
  fh=open(per_chain_repertoire_statistics_file, "w")
  fh.write(out)
  fh.close()
  return()

def Get_locations(ref_locations):
  fh=open(ref_locations,"r")
  locations = {}
  for l in fh:
    l=l.strip().split()
    locations[l[0]] = l[1]
  fh.close()
  return(locations)

###########################
dir = sys.argv[1]             ### output directory
id = sys.argv[2]              ### identifier for output files
seq_file= sys.argv[3]         ### input fasta
protein_file = sys.argv[4]    ### protein fasta (does not matter if read number index is not correct)
gene = sys.argv[5]            ### either IGH, IGL, IGK, TRA, TRB, TRD, TRG
species = sys.argv[6]         ### HOMO_SAPIENS...
cluster_file = sys.argv[7]    ### cluster file (output from Read_processing_and_quality.py program)
command = sys.argv[8]         ### Comma separated list: ANNOTATE = Generate raw V-J annotation files, STATISTICS = Get population statistics
forward_primer_group = sys.argv[9] ### Forward primers used (basically to disinguish between ISO_DD and IsoTyper)

########################### Files for annotation generation
annot_file = dir+"TMP/Annotation_"+id+".txt"
tmp_file =dir+"TMP/Tmp_annotation_"+id
network_statistics = dir+"Network_statistics_"+id+".txt"
cluster_statistics = dir+"Cluster_statistics_"+id+".txt"
gene_freq_file = dir+"Gene_frequencies_"+id+".txt"
cluster_properties_file = dir+"TMP/Cluster_properties_file_"+id+".txt"
constant_region_count_file = dir+"Constant_region_counts_"+id+".txt"
cluster_size_distribution = dir+"Distribution_cluster_sizes_"+id+".txt"
vertex_size_distribution = dir+"Distribution_vertex_sizes_"+id+".txt"
per_chain_repertoire_statistics_file = dir+"IsoTyper_chain_repertoire_statistics_file_"+id+".txt"
########################## Reference files
loc = "/nfs/users/nfs_r/rbr1/IMMUNOLOGY/NETWORK_METHODS/STANDARD_BCR_PROCESSING/"
loc = "/users/immune-rep/mfj169/REPERTOIRE_ANALYSES/"
refv = loc+"LIBRARY/Reference_nn_"+species+"_"+gene+"V.fasta"
refd = loc+"LIBRARY/Reference_nn_"+species+"_"+gene+"D.fasta"
refj = loc+"LIBRARY/Reference_nn_"+species+"_"+gene+"J.fasta"
refvp= loc+"LIBRARY/Reference_protein_"+species+"_"+gene+"V.fasta"
refjp= loc+"LIBRARY/Reference_protein_"+species+"_"+gene+"J.fasta"
CDR3_end_file = loc+"LIBRARY/CDR3_end_"+species+"_imgt.txt"
codon_file = loc+"LIBRARY/Codon_table2.txt"
########################## Get commands
pwd = "/well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/"#commands.getoutput("pwd")
ref_locations = pwd+"/Locations_of_called_programmes.txt"
locations = Get_locations(ref_locations)
command = command.split(",")
#constant_region = "TRUE"
constant_region = "FALSE"
isotyper_primer_set = "INNER"
print constant_region, isotyper_primer_set

######################### Commands
Intialise_files(dir)
if("ANNOTATE" in command):
  a = 1
  if(species in ["HOMO_SAPIENS","MACACA_MULATTA","MUS_MUSCULUS","LLAMA_GLAMA"]):
    Local_immune_repertoire_annotator = locations["Local_immune_repertoire_annotator"]
    command1 = "python "+Local_immune_repertoire_annotator+" "+dir+"TMP/ "+id+" "+seq_file+" "+gene+" "+species
    commands.getoutput(command1)
    Get_gene_frequencies(annot_file, gene_freq_file,gene,id)
    Get_constant_region_distribution(seq_file,constant_region_count_file,id,annot_file)
if("STATISTICS" in command):
  Initalise_file(network_statistics)
  Initalise_file(cluster_statistics)
  Get_network_statistics(cluster_file, id, dir,network_statistics,gene,species,cluster_size_distribution,vertex_size_distribution)
  Get_cluster_stats(annot_file, cluster_file, id, dir,cluster_statistics,gene,species,seq_file,loc)
  Get_network_statistics_per_chain(cluster_file, id, dir,per_chain_repertoire_statistics_file,isotyper_primer_set)

