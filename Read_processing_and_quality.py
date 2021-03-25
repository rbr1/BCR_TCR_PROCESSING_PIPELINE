# ImmuneReceptor_PROCESSING_PIPELINE developed by Rachael Bashford-Rogers (2020)
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
import commands
import sys
from operator import itemgetter
import networkx as nx
import re

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

def fuzzy_substring(needle, haystack):
    """Calculates the fuzzy match of needle in haystack,
    using a modified version of the Levenshtein distance
    algorithm.
    The function is modified from the levenshtein function
    in the bktree module by Adam Hupp"""
    m, n = len(needle), len(haystack)
    # base cases
    if m == 1:
        return not needle in haystack
    if not n:
        return m
    row1 = [0] * (n+1)
    for i in range(0,m):
        row2 = [i+1]
        for j in range(0,n):
            cost = ( needle[i] != haystack[j] )

            row2.append( min(row1[j+1]+1, # deletion
                               row2[j]+1, #insertion
                               row1[j]+cost) #substitution
                           )
        row1 = row2
    return min(row1)

def Intialise_files(dir):
  dirs_to_add=["FASTQ_FILES/", "ORIENTATED_SEQUENCES/", "ORIENTATED_SEQUENCES/TMP/", "ORIENTATED_SEQUENCES/NETWORKS/"]
  for d in dirs_to_add:
    c=commands.getoutput("ls "+dir+d)
    if(c.count("No such file or directory")==1):commands.getoutput("mkdir "+dir+d)
  return()

def QC_samples(dir,gene,id,source,length,species,barcode_group,quasr_location,readjoin_location):
  pre = source
  reads1 = source+"R1_001.fastq"
  reads2 = source+"R2_001.fastq"
  #reads1 = source
  #reads2 = source.replace("1.fq","2.fq")
  reads1 = source+"1.fastq"
  reads2 = source+"2.fastq"
  #reads1 = source.replace("XXXXX","1")
  #reads2 = source.replace("XXXXX","2")
 # os.system("gunzip "+reads1+".gz")
 # os.system("gunzip "+reads2+".gz")
  if(1==1):
    print reads1
    print reads2
    if(1==1):
      command0 = readjoin_location+" "+reads1 +" "+reads2 +" -O -o Joined_"+id+" -d "+dir+"FASTQ_FILES/  --min-overlap=15  --max-overlap=280"
      os.system(command0)
      threshold,length = "30","100"
      command1 = "java -jar "+quasr_location+"qualityControl.jar -f "+dir+"FASTQ_FILES/Joined_"+id+".extendedFrags.fastq -o "+dir+"FASTQ_FILES/Sequences_"+id+"_1 -m "+threshold+" -l "+length
      command2 = "cat "+dir+"FASTQ_FILES/Sequences_"+id+"_1.qc.fq | perl -e '$i=0;while(<>){if(/^\@/&&$i==0){s/^\@/\>/;print;}elsif($i==1){s/\./N/g;print;$i=-3}$i++;}' > "+dir+"FASTQ_FILES/Sequences_"+id+"_1.fasta"
      print command1
      os.system(command1)
      os.system(command2)
  return()

def Get_paired_reads_overlapping(file1, file2, outfile,gene,paired,id,method):
  (rc)= Init_rc()
  (seqs1)=Get_sequences(file1)
  (seqs2)=Get_sequences(file2)
  print "Forward reads:", len(seqs1),"Reverse reads:", len(seqs2)
  out,ind = '',0
  fh = open (outfile, "w")
  fh.close()
  tot,f_joining,total,length,fail_no_pair=0,0,0,20,0
  if(method in ["5RACE"]):length = 20
  if(gene=="LAMBDA" or gene=="KAPPA" or gene =="IGK"):length=10
  file_prefix = dir+"FASTQ_FILES/Seqs_barcode_split_"+barcode_group+"_"
  print  length,"length overlapping"
  for id in seqs1:
    if (id in seqs2):
      total=total+1
      if(min([len(seqs1[id]), len(seqs2[id])])>250):
        (seq,failed)=Join_reads(seqs1[id], seqs2[id],rc,length)
        print id, seq, failed
        if (failed!=0 or len(seq) < 300):
          (seq,failed)=Join_reads(seqs1[id], seqs2[id],rc,6)
          if (failed!=0 or len(seq) < 300):
            ind,tot = ind+1,tot+1
            out = out+">"+id+"\n"+seq+"\n"
        if (failed==0 and len(seq)>300):
          seq2 = Reverse_comp(seqs2[id], rc)
          (seq,failed)=Join_reads(seqs1[id], seq2,rc,length)
          if (failed==0 and len(seq)>180):
            ind,tot = ind+1,tot+1
            out = out+">"+id+"\n"+seq+"\n"
          else:
            f_joining=f_joining+1
            #print seqs1[id], seqs2[id]
        if (ind >100):
          Write_out(out, outfile)
          out,ind = '',0
    else:
      fail_no_pair=fail_no_pair+1
  Write_out(out, outfile)
  del out,ind
  print id+"\tTotal sequences:",total,"\tNumber of Failed sequences:",f_joining, "\tPercentage sequences failed:", f_joining*100.0/total,"%\tTotal remaining:",tot, "\tNo pairs:",fail_no_pair
  return()

def Get_barcodes(barcode, barcode_file,rc):
  barcode = barcode.split(",")
  fh=open(barcode_file,"r")
  barcodes,allbarcodes,barcode_stem,barcode_only,barcode_only_rc = [],[],[],[],[]
  inverse,barcode_seq = {},{}
  rc_barcodes = []
  for l in fh:
    if(l[0]!="#"):
      l=l.strip().split()
      allbarcodes.append(l[1])
      if(l[0]!="STEM"):
        barcode_seq[l[0]]=l[1]
        barcodes.append(l[1])
        rc_barcodes.append(Reverse_comp(l[1], rc))
        barcode_only.append(l[2])
        barcode_only_rc.append(Reverse_comp(l[2], rc))
        sq = [l[1],Reverse_comp(l[1], rc),l[2],Reverse_comp(l[2], rc)]
        for s in sq: inverse[s] = l[0]
      elif(l[0]=="STEM"):barcode_stem.append([l[1],l[1][0:len(l[1])/2], l[1][len(l[1])/2:len(l[1])]])
  fh.close()
  return(barcodes,allbarcodes,barcode_stem, rc_barcodes,barcode_only,barcode_only_rc,inverse,barcode_seq)

def Init_rc():
  b1 = ["A","T","G","C","N","."]
  b2 = ["T","A","C","G","N","."]
  rc = {}
  for i in range(0,6):
    rc[b1[i]]=b2[i]
  return(rc)

def Reverse_comp(seq, rc):
  s = ''
  l=len(seq)
  for i in range(0,l):
    j=l-i-1
    if(seq[j] not in rc):
      print seq
    else:
      s=s+rc[seq[j]]
  return(s)

def Convert_sequence_headers(Seq_file1, Tmp_file):
  fh = open(Tmp_file,"w")
  fh.close()
  fh = open (Seq_file1, "r")
  out, ind = '',0
  for header,sequence in fasta_iterator(fh):
    header= header.replace(":","").split()[0].replace("-","")
    header = header.split("#")[0].split("/")[0]
    out=out+">"+header+"\n"+sequence+"\n"
    ind = ind+1
    if(ind>500):
      Write_out(out, Tmp_file)
      out, ind = '',0
  fh.close()
  Write_out(out, Tmp_file)
  out, ind = '',0
  return()

def Get_sequences(file):
  command = "gunzip "+file
  commands.getoutput(command)
  fh = open (file, "r")
  seqs = {}
  for header,sequence in fasta_iterator(fh):
    header = header.replace(":","").split()[0].replace("-","")
    header = header.split("#")[0].split("/")[0]
    seqs[header]=sequence
    #if(len(seqs)>20000):break
  fh.close()
  return(seqs)

def Split_reads_and_join(file1, file2, outfile,gene,paired,id,method,primer_file,barcode_group,other,dir,ref_const):
  print method,paired
  if(method=="Multiplex_FORWARD_BARCODE_GROUPED" or method=="Multiplex_REVERSE_BARCODE_GROUPED"):
    Split_reads_and_join_pre(file1, file2, outfile,gene,paired,id,method,primer_file,barcode_group,other,dir,ref_const)
  elif(paired=="OVERLAPPING" and method in ["Multiplex"]):
    #Get_paired_reads_overlapping(Seq_file1, Seq_file2, Tmp_file,gene,paired,id,method)
    print Tmp_file
    print Seq_file1
    Convert_sequence_headers(Seq_file1, Tmp_file)
  return()

def Check_split(outfile, barcode_group,other):
  rc = Init_rc()
  barcode, barcode_file = other.split(",")[0], other.split(",")[1]
  barcodes,allbarcodes,barcode_stem,rc_barcodes,barcode_only,barcode_only_rc = Get_barcodes(barcode, barcode_file,rc)
  allbarcodes_rc = []
  for b in allbarcodes:
    allbarcodes_rc.append(Reverse_comp(b,rc))
  out,ind = '',0
  outfile1 = outfile.replace(".fasta","_CHECKED.fasta")
  fh=open(outfile1,"w")
  fh.close()
  fh = open(outfile,"r")
  seqs = {}
  for header,sequence in fasta_iterator(fh):
    seqs[header] = sequence
  fh.close()
  for header in seqs:
    sequence = seqs[header]
    score = []
    for i in range(0,len(allbarcodes)):
      p1 = len(Get_match(allbarcodes[i], sequence))
      p2 = len(Get_match(allbarcodes_rc[i], sequence))
      score.append(p1+p2)
    if(sum(score)<=1):
      out=out+">"+header+"\n"+sequence+"\n"
      ind = ind+1
      if(ind>500):
        Write_out(out, outfile1)
        out, ind = '',0
  fh.close()
  Write_out(out, outfile1)
  return()

def Join_reads(s1, s2,rc,length):
  seq=''
  (s2)=Reverse_comp(s2,rc)
  (l1,l2)=(len(s1),len(s2))
  failed = 1
  for i in range(0,100):
    ind=(i*5)+5
    if(ind<(l1-length)):
      (s11,s22,p,overlap)=Trim(s1,s2,l1,l2,ind,length)
      if(p==1):
        seq = s1[0:s1.index(overlap)]+s2[(s2.index(overlap)):l2].lower()
        if(len(seq)>200):
          failed = 0
          break
  if(failed==1):
    seq = s1+"..."+s2
  return(seq,failed)

def Trim(s1,s2,l1,l2,indent,length):
  (i,p)=(indent,0)
  sample1=s1[i:i+length]
  index=s2.find(sample1)
  s1a, s2a = s1, s2
  if (index!=-1):
    if(index>i):
      s2a=s2a[index-i:l2]
    else:
      s1a=s1a[i-index:l1]
    min_len=min([len(s1),len(s2)])
    s1a=s1a[0:min_len]
    s2a=s2a[0:min_len]
    p=1
  return (s1a, s2a, p,sample1)

def Get_match(primer, seq):
  loc = []
  if (seq.count(primer)!=0):
    for m in re.finditer(primer, seq):
      loc = loc+[m.start()]
  return(loc)

def  Get_gene_ref(ref_const,kmer):
  fh=open(ref_const, "r")
  ref_lib = Tree()
  for header,sequence in fasta_iterator(fh):
    sequence = sequence.upper()
    header = header.split("*")[0]
    for i in range(kmer, len(sequence)):
      k = sequence[i-kmer: i]
      ref_lib[k][header].value = 1
  fh.close()
  return(ref_lib)

def Check_sequences_for_genes(s, ref_lib,kmer):
  scores = {}
  for i in range(kmer, len(s)):
    k =s[i-kmer: i]
    if(k in ref_lib):
      for const in ref_lib[k]:
        if(const in scores):scores[const] = scores[const]+1
        else:scores[const]=1
  maxi, maxi_iso = 0,''
  for c in scores: 
    if(scores[c]>maxi):  maxi, maxi_iso =scores[c], c
  return(maxi, maxi_iso)

def Split_reads_and_join_pre(file1, file2, outfile,gene,paired,sample,method,primer_file,barcode_group,other,dir,ref_const):
  kmer = 12
  direction = 1
  if(method =="Multiplex_FORWARD_BARCODE_GROUPED"):direction = 1
  if(method =="Multiplex_REVERSE_BARCODE_GROUPED"):direction = 2
  ref_lib = Get_gene_ref(ref_const,kmer)
  fh=open(outfile,"w")
  fh.close()
  rc = Init_rc()
  print "\t",other
  barcode, barcode_file = other.split(",")[0], other.split(",")[1]
  barcodes,allbarcodes,barcode_stem,rc_barcodes,barcode_only,barcode_only_rc,inverse,barcode_seq = Get_barcodes(barcode, barcode_file,rc)
  (rc)= Init_rc()
  (seqs1)=Get_sequences(file1)
  print len(seqs1)
  out,ind = '',0
  fh = open (outfile, "w")
  fh.close()
  tot,joining,total,length,primer_match=0,0,0,8,0
  if(gene=="LAMBDA" or gene=="KAPPA" or gene =="IGK"):length=30
  count_joined, count_barcode, count_barcode_saved = 0,0,0
  file_prefix = dir+"FASTQ_FILES/Seqs_barcode_split_"+sample+"_"
  for bc in barcodes:
    outfile = file_prefix+str(inverse[bc])+".fasta"
    fh=open(outfile,"w")
    fh.close()
  seqs_out = {}
  for id in seqs1:
    total=total+1
    seq = seqs1[id]
    if(len(seq)>250):
      score1, gene1 = Check_sequences_for_genes(seq, ref_lib,kmer)
      score2, gene2 = Check_sequences_for_genes(Reverse_comp(seq,rc), ref_lib,kmer)
      if(max([score1, score2])>8): ### 40 for the longer IsoTyper sequences
        barcode_matches = {}
        for primer in barcodes:
          p = Get_match(primer, seq)
          if(len(p)==0):
            seq = Reverse_comp(seq,rc)
            p = Get_match(primer, seq)
            if(len(p)!=0):
              if(p[0]<len(seq)/2):
                bc = inverse[primer]
                barcode_matches[primer+"|"+bc+"|"+str(p[0])]=1
          if(len(p)!=0):
            if(p[0]<len(seq)/2):
              bc = inverse[primer]
              barcode_matches[primer+"|"+bc+"|"+str(p[0])]=1
        if(len(barcode_matches)>0):
          if(len(barcode_matches)==1):
            for b in barcode_matches: break
            bc = b.split("|")[1]
            count_barcode = count_barcode+1
            if(seq.count(barcode_seq[bc])==0):
              if(Reverse_comp(seq,rc).count(barcode_seq[bc])==1):
                seq = Reverse_comp(seq,rc)
            if(direction==2):seq = Reverse_comp(seq,rc)
            out1=">"+id+"\n"+seq+"\n"
            if(bc in seqs_out): seqs_out[bc] = seqs_out[bc]+out1
            else:seqs_out[bc]=out1
            ind,primer_match = ind+1,primer_match+1
            if(seqs_out[bc].count(">")>200):
              outfile = file_prefix+str(bc)+".fasta"
              Write_out(seqs_out[bc], outfile)
              seqs_out[bc] = ''
        else:
          for prim in barcode_stem:
            p = Get_match(prim[0], seq)
            if(len(p)==0):
              seq = Reverse_comp(seq,rc)
              p = Get_match(prim[0], seq)
            if(len(p)!=0):
              start = p[0]+len(prim[0])-3
              start = p[0]+len(prim[0])-6
              for i in range(0,len(barcodes)):
                bc_halves = [barcodes[i][2:len(barcodes[i])/2], barcodes[i][len(barcodes[i])/2:len(barcodes[i])-2]]
                bc_shift = [2,len(barcodes[i])/2]
                test = seq[start: start+len(barcodes[i])+10]
                for j in range(0,len(bc_halves)):
                  if(test.count(bc_halves[j])!=0):
                    start = barcodes[i].index(barcode_only[i])
                    p=test.index(bc_halves[j])-bc_shift[j]
                    sub_test = test[p+start:p+len(barcode_only[i])+start]
                    a = fuzzy_substring(sub_test,barcode_only[i])
                    if(a<=1):
                      primer = barcodes[i]
                      bc = inverse[primer]
                      barcode_matches[primer+"|"+bc+"|"+str(0)]=1
                    else:
                      seq1 = Reverse_comp(seq,rc)
                      a = fuzzy_substring(sub_test,barcode_only[i])
                      if(a<=1):
                        primer = barcodes[i]
                        bc = inverse[primer]
                        barcode_matches[primer+"|"+bc+"|"+str(0)]=1
                        seq =seq1
          if(len(barcode_matches)>0):
            if(len(barcode_matches)==1):
              for b in barcode_matches: break
              bc = b.split("|")[1]
              count_barcode = count_barcode+1
              if(direction==2):seq = Reverse_comp(seq,rc)
              out1=">"+id+"\n"+seq+"\n"
              if(bc in seqs_out): seqs_out[bc] = seqs_out[bc]+out1
              else:seqs_out[bc]=out1
              ind,primer_match = ind+1,primer_match+1
              if(seqs_out[bc].count(">")>200):
                outfile = file_prefix+str(bc)+".fasta"
                Write_out(seqs_out[bc], outfile)
                seqs_out[bc] = ''
  for bc in seqs_out:
    outfile = file_prefix+str(bc)+".fasta"
    Write_out(seqs_out[bc], outfile)
    seqs_out[bc] = ''
  del seqs_out
  print id+"\tTotal sequences:",total,"\tNumber of joined sequences:",joining, "Number of barcodes matched:",primer_match,"% retained:",primer_match*100.0/total, "IDs saved",count_barcode_saved
  return()

def Single_J_barcoded_trimming_clustered(forward, reverse,  barcoded_j, barcoded_v, rc,Tmp_file,Fail_file,Output_trim,primer_tag_file,tmp_file,gene,paired,species,primer_file,primer_tag_file_count,threshold_BARCODE,ref_const,inside,v_ref,cd_hit_directory):
  Read_untrimmed_file_single(rc,Tmp_file,Fail_file,Output_trim, gene,paired,species,primer_file,primer_tag_file,tmp_file,primer_tag_file_count,forward, reverse,v_ref)
  Check_barcodes_MALBAC(primer_tag_file_count,primer_tag_file,Fail_file,Output_trim,threshold_BARCODE,cd_hit_directory,tmp_file)
  Separate_sequences(primer_tag_file_count,primer_tag_file, Output_trim,ref_const)
  return()

def Assess_gene_score(consensus, word_dict,threshold):
  scores,max_score = [],0
  for hk in word_dict:
    w = word_dict[hk]
    score = [i for i in range(len(w)) if consensus.count(w[i][0])!=0]
    if(len(score)>threshold/2):
      if(len(score)>=max_score):
        consensus1 = consensus
        if(hk.count("HOUSEKEEPING")==0 and len(score)>threshold/2):
          start = [consensus.index(w[score[i]][0])-w[score[i]][1] for i in range(len(score))]
          start = max(set(start), key=start.count)
          consensus1 = consensus[0:start]
        if(max_score==len(score)):
          scores.append([len(score),hk,consensus1])
        else:
          scores = []
          scores.append([len(score),hk,consensus1])
        max_score = len(score)
  return(scores)

def Separate_sequences(primer_tag_file_count,primer_tag_file, Output_trim,ref_const):
  fh=open(ref_const, "r")
  word_length = 8
  word_dict = {}
  for header,sequence in fasta_iterator(fh):
    header,sequence= header.replace("|",",").split("*")[0],sequence.upper()
    words= []
    for i in range(word_length,len(sequence)):
      words.append([sequence[i-word_length:i],i-word_length])
    word_dict[header] = words
  fh.close()
  region_primer_type = {}
  fh=open(primer_tag_file_count,"r")
  for l in fh:
    if(l[0]!="#"):
      l=l.strip().split()
      const_reg,seq_type,rev_primer = l[len(l)-1],l[4],l[6]
      region_primer_type[l[0].split(":")[0]] = [const_reg,seq_type,rev_primer]
  fh.close()
  fh=open(primer_tag_file,"r")
  hk_score,non_hkscore = [],[]
  seqs,classifications = Tree(),{}
  for l in fh:
    if(l[0]!="#"):
      l=l.strip().split()
      if(l[6]=="YES"):
        consensus = l[7]
        id = l[0].split("__")[0]
        seqs[consensus][id].value = 1
  fh.close()
  print len(seqs)
  out,ind = '',0
  seq_info, seq_uniq = {},Tree()
  classes = {}
  threshold = 5
  for s in seqs:
    for id in seqs[s]:
      break
    scores= Assess_gene_score(s,word_dict,threshold)
    if(len(scores)>0):
      if(scores[0][0]>threshold): ## for non-housekeeping genes, trimmed sequences to remove constant region
        if(len(scores)==1):
          if(len(scores[0][2])>50):seq = scores[0][2]
          else:seq = s
          seq_uniq[s][id].value = 1
          seq_info[id] = [len(seqs[s]), seq,scores[0][1]]
          classes[scores[0][1]] = 1
  del seqs
  print len(classes)
  classes_all = []
  for c in classes:
    classes_all.append(c)
  classes_all.sort()
  classes_order ={}
  for i in range(0,len(classes_all)):
    classes_order[classes_all[i]] = i
  header = "_".join(classes_all)
  for s in seq_uniq:
    if(len(seq_uniq[s])>0):
      f=[0]*len(classes_all)
      for id in seq_uniq[s]:
        c = seq_info[id][2]
        f[classes_order[c]] = f[classes_order[c]]+seq_info[id][0]
      out=out+">"+id.split("__")[0]+"__"+"_".join(map(str, f))+"|"+header+"\n"+s+"\n"
      ind = ind+1
      if(ind>100):
        Write_out(out, Output_trim)
        out, ind = '',0
  Write_out(out, Output_trim)
  return()

def Check_barcodes_MALBAC(primer_tag_file_count,primer_tag_file,Fail_file,Output_trim,threshold,cd_hit_directory,tmp_file):
  outs = ['','#ID\tnumber_of_reads\ttotal_reads_with_BC\tJ_barcode\tV_barcode\tbp_mismatches_to_consensus\tBC_accepted\tconsensus\tsequence\n']
  files = [Output_trim,primer_tag_file]
  for i in range(0,len(files)):
    fh=open(files[i],"w")
    fh.write(outs[i])
    fh.close()
  fh=open(primer_tag_file_count,"r")
  seqs = Tree()
  for l in fh:
    if (l[0]!="#"):
      l=l.strip().split()
      v_tag, j_tag, sequence, header = l[2],l[1],l[3],l[0]
      seqs[j_tag+"\t"+v_tag][sequence][header].value = 1
  fh.close()
  print len(seqs), "Unique tags"
  min_depth = (1.0/(1-threshold))
  total_tags,ind,passed_seqs_total = 0,0,0
  fail_less_than_threshold = 0
  outp=''
  for t1 in seqs:
    t,u_seq, u_freq,u_header = 0,[],[],[]
    for s in seqs[t1]:
      f = 0 
      for h in seqs[t1][s]:
        f = f+int(h.split(":")[1])
      total_tags,u_seq, u_freq, u_header, ind =total_tags+f, u_seq+[s], u_freq+[f], u_header+[h], ind+1
    f = sum(u_freq)
    h = h.split(":")[0].split("__")[0]+"__"+str(sum(u_freq))
    if(sum(u_freq)>20):print "BC group size:\t",len(u_freq),t1.split() 
    if(len(u_freq)==1): 
      passed_seqs_total = passed_seqs_total+f
      outp=outp+h+"\t"+str(f)+"\t"+str(f)+"\t"+t1+"\t0\tYES\t"+s+"\t"+s+"\n"
    elif(len(u_freq)<500):
      #print "BC group size:\t",len(u_freq),t1.split(), u_freq#, u_seq
      if(max(u_freq)>sum(u_freq)*threshold):
        nz = [i for i in range(len(u_freq)) if u_freq[i]!=max(u_freq)]
        passed_seqs_total = passed_seqs_total+f
        consensus = u_seq[nz[0]]
        outp = outp+"\t".join(map(str, [h, len(u_freq), sum(u_freq),t1,0,"YES",consensus, consensus ]))+"\n"
      elif(len(u_freq)>15 and len(u_freq)<200): ## clustering first then alignment
        out_cluster,ids = '',{}
        for i in range(0,len(u_seq)):
          out_cluster = out_cluster+">"+u_header[i]+"\n"+u_seq[i].replace("...","")+"\n"
          ids[u_header[i]] = [u_seq[i], u_freq[i]]
        consensus,pass_consensus = Get_consensus_sequence_large(out_cluster, Fail_file,len(u_seq[i]), ids,sum(u_freq),threshold,tmp_file,cd_hit_directory)
        if(consensus!='' and pass_consensus==1):
          outp = outp+"\t".join(map(str, [h, len(u_freq), sum(u_freq),t1,0,"YES",consensus, consensus ]))+"\n"
          passed_seqs_total = passed_seqs_total+f
        else:fail_less_than_threshold = fail_less_than_threshold+1
      else: 
        consensus, pass_consensus=Get_consensus_sequence_cluster(u_seq, u_freq,tmp_file,threshold)
        if(consensus.count("_")==0 and pass_consensus==1):
          outp = outp+"\t".join(map(str, [h, len(u_freq), sum(u_freq),t1,0,"YES",consensus, consensus ]))+"\n"
          passed_seqs_total = passed_seqs_total+f
        else:
          fail_less_than_threshold = fail_less_than_threshold+1
    if(ind>200):
       Write_out(outp, primer_tag_file)
       outp, ind = '',0
  Write_out(outp, primer_tag_file)
  outp, ind = '',0
  print total_tags, passed_seqs_total, fail_less_than_threshold
  return()

def Get_consensus_sequence_cluster(u_seq, u_freq,tmp_file,threshold):
  out=''
  for i in range(0,len(u_seq)):
    out=out+">"+str(i)+"\n"+u_seq[i].replace("...","")+"\n"
  fh=open(tmp_file+"txt", "w")
  fh.write(out)
  fh.close()
  insert = ''
  if( len(u_seq) > 2000):insert = '--parttree'
  mafft = "/apps/well/mafft/7.149/bin/mafft "
  command1 = mafft+" --retree 2 "+insert+" "+tmp_file+"txt > "+tmp_file+"aligned"
  commands.getoutput(command1)
  fh=open(tmp_file+"aligned","r")
  max_seqs={}
  for header,sequence in fasta_iterator(fh):
    max_seqs[sequence.upper()]=int(header)
  fh.close()
  bases=["A","T","G","C","-"]
  base_dict = {}
  for b in range(0,len(bases)):
    base_dict[bases[b]]=b
  consensus = ''
  start,pass_consensus = 0,1
  for i in range(0,len(sequence)):
    f = [0]*len(bases)
    for s in max_seqs:
      if(s[i] not in base_dict):
        print i, s[i], max_seqs
      f[base_dict[s[i]]]=f[base_dict[s[i]]]+u_freq[max_seqs[s]]
    if(f[4]==0):start==1
    if(max(f)*1.0/sum(f) >= threshold):
      if(bases[f.index(max(f))]!="-"):
        consensus=consensus+bases[f.index(max(f))]
    else:
      pass_consensus = 0
      if(start==1):
        for j in range(0,5):
          if(f[j]!=0):
            consensus=consensus+"|"+bases[j]+":"+str('%s' % float('%.3g' % (f[j]*1.0/sum(f))))
        consensus=consensus+"_"
      else:
        f = f[0:4]
        if(bases[f.index(max(f))]!="-"):
          consensus=consensus+bases[f.index(max(f))]
        else:
          for j in range(0,4):
            if(f[j]!=0):
              consensus=consensus+"|"+bases[j]+":"+str('%s' % float('%.3g' % (f[j]*1.0/sum(f))))
            consensus=consensus+"_"
  return(consensus,pass_consensus)

def Get_consensus_sequence_large(out_cluster, Fail_file,l1,ids,sum_u_freq,threshold,tmp_file,cd_hit_directory):
  fh=open(Fail_file,"w")
  fh.write(out_cluster)
  fh.close()
  Cluster_i(Fail_file,Fail_file+"cls",(l1-3.0)/l1,cd_hit_directory)
  cluster = Tree()
  fh=open(Fail_file+"cls.bak.clstr","r")
  for l in fh:
    l=l.strip().split()
    clust, id = l[0],l[2].replace("...","").replace(">","")
    cluster[clust][id].value = 1
  fh.close()
  max_s, max_clust = 0,''
  for c in cluster:
    f = 0
    for id in cluster[c]:
      f = f+ids[id][1]
    if(max_s<f):max_clust,max_s = c,f
  consensus,pass_consensus = '',0
  if(sum_u_freq*threshold<max_s):
    seq_align,s_freq = [],[]
    for id in cluster[max_clust]:
      seq_align.append(ids[id][0])
      s_freq.append(ids[id][1])
    consensus = Get_consensus_sequence(seq_align,s_freq,tmp_file,threshold)
    if(consensus.count("_")==0 and len(consensus)>3):pass_consensus=1
    else:pass_consensus = 0
  return(consensus, pass_consensus)

def Read_untrimmed_file_single(rc,Tmp_file,Fail_file,Output_trim, gene,paired,species,primer_file,primer_tag_file,tmp_file,primer_tag_file_count,forward, reverse,v_ref):
  for f in [primer_tag_file_count]:
    fh=open(f, "w")
    fh.close()
  fh = open (Tmp_file, "r")
  seqs = Tree() 
  minl,maxl=110,1000000
  if(gene=="HEAVY" or gene == "IGH"):minl= 120
  if(gene=="KAPPA" or gene=="IGK"):minl= 110 
  if(gene=="LAMBDA" or gene == "IGL"):(minl, maxl)= (90, 150)
  J_found,v_found,tot,indexing=0,0,0,0
  seqs1,t=Tree(),0
  for header,seq in fasta_iterator(fh):
    seq=seq.upper()
    if(seq.count("...")==0 and len(seq)>180):
      seqs1[seq][header].value=1
      t=t+1
  out,ind = "#ID\tJ_tag\tV_tag\tSequence\ttype_rev\ttype_for\tprimer_rev\tprimer_for\n", 0
  total,pass_r, pass_f, pass_all = 0,0,0,0
  print "start"
  for seq in seqs1:
    for header in seqs1[seq]:
      break
    header = header+":"+str(len(seqs1[seq]))
    number = len(seqs1[seq])
    passes,j_tag,v_tag=0,'',''
    type_rev, type_for = '',''
    primer_rev,primer_for = '',''
    total = total+1
    #if(seq.count("CAGCACGCTTCAGGCT")!=0): print header
    for i in range(0,len(reverse)):
      pj=Get_match(reverse[i][2], seq)
      if(max(pj+[-1])==-1 or len(pj)==0):
        seq = Reverse_comp(seq, rc) 
        pj=Get_match(reverse[i][2], seq)
      if(max(pj+[-1])!=-1):
        bc_len = len(reverse[i][3])#len(reverse[i][2])/5
        if(pj[0]>bc_len-3):
          #j_tag = seq[pj[0]-bc_len:pj[0]]
          j_tag = seq[pj[0]+len(reverse[i][2]): pj[0]+len(reverse[i][2])+bc_len]
          if(len(j_tag)>bc_len/2):
            seq = seq[pj[0]+bc_len:len(seq)]
            seq = Reverse_comp(seq, rc)
            type_rev,primer_rev = reverse[i][1], reverse[i][5]
            passes = 1
            break
    if(passes!=1):
      for i in range(0,len(reverse)):
        words =reverse[i][6]
        p=[]
        for w in words:
          pj=Get_match(w[0], seq)
          if(max(pj+[-1])!=-1 and len(pj)>0):
            if(pj[0]<len(seq)/2):p=p+[max([0,pj[0]-w[1]])]
            else:pj[0]=-1
        if(len(p)>1):
          pj = max(set(p), key=p.count)
          bc_len = len(reverse[i][3])
          if(pj>min([bc_len-3])):
            j_tag = seq[pj-bc_len:pj+1]
            if(len(j_tag)>bc_len/2):
              seq = seq[pj+bc_len:len(seq)]
              seq = Reverse_comp(seq, rc)
              type_rev,primer_rev = reverse[i][1], reverse[i][5]
              passes = 1
              break
      #if(passes!=1):
      #  seq = Reverse_comp(seq, rc) 
      #  for i in range(0,len(reverse)):
      #    words =reverse[i][6]
      #    p=[]
      #    for w in words:
      #      pj=Get_match(w[0], seq)
      #      if(max(pj+[-1])!=-1 and len(pj)>0):
      #        if(pj[0]<len(seq)/2):p=p+[max([0,pj[0]-w[1]])]
      #        else:pj[0]=-1
      #    if(len(p)>1):
      #      bc_len = len(reverse[i][3])
      #      pj = max(set(p), key=p.count)
      #      if(pj>min([bc_len-3])):
      #        j_tag = seq[pj-bc_len:pj+1]
      #        print header, j_tag, seq[pj-bc_len-10: pj+len(reverse[i][2])+bc_len]
      #        if(len(j_tag)>bc_len/2):
      #          seq = seq[pj+bc_len:len(seq)]
      #          seq = Reverse_comp(seq, rc)
      #          type_rev,primer_rev = reverse[i][1], reverse[i][5]
      #          passes = 1
      #          break
    if(passes==1):
      pass_r = pass_r+1
      passes = 2
      v_tag = '-'
      type_for, primer_for = '-','-'
      #for i in range(0,len(forward)):
      #  pv=Get_match(forward[i][0], seq)
      #  if(max(pv+[-1])!=-1):
      #    v_tag='-'
      #    seq = seq[pv[0]+len(forward[i][0])-1 :len(seq)]
      #    type_for, primer_for = forward[i][1],forward[i][2]
      #    passes = 2
      #    break
      #if(passes!=2):
      #  for i in range(0,len(forward)):
      #    words =forward[i][3]
      #    p=0 
      #    for w in words:
      #      pv=Get_match(w[0], seq) 
      #      if(max(pv+[-1])!=-1 and len(pv)>0):
      #        if(pv[0]<len(seq)/2):p=p+1
      #        else:pv[0]=-1
      #    if(max(pv+[-1])!=-1 and len(pv)>0 and p>1):
      #      v_tag='-'
      #      seq = seq[pv[0]+len(forward[i][2])-1 :len(seq)]
      #      type_for, primer_for = forward[i][1],forward[i][2]
      #      passes = 2
      #      break
      #if (passes!=2):
      #  p=[]
      #  for i in range(0,len(v_ref)):
     #     words= v_ref[i][3]
     #     p1 = []
     #     for w in words:
     #       pv=Get_match(w[0], seq)
     #       if(max(pv+[-1])!=-1 and len(pv)>0):
     #         if(pv[0]<len(seq)/2):p=p+[pv[0]]
     #     if(len(p1)>0):p1 = p1+[min(p1)]
     #   if(len(p)>=2):
     #     pv = min(p)
     #     v_tag='-'
     #     seq = seq[pv:len(seq)]
     #     type_for, primer_for = v_ref[i][1],v_ref[i][2]
     #     passes = 2
    if (passes ==2):
      pass_f = pass_f+1
      if(len(seq)>minl and len(seq)<maxl):
        if(seq.count("N")==0):
          #if(type_rev==type_for):
          out=out+"\t".join(map(str,[header, j_tag, v_tag, seq, type_rev,type_for,primer_rev,primer_for]))+"\n"
          pass_all = pass_all+1
          ind = ind+1
          if(ind>200):
            Write_out(out, primer_tag_file_count)
            out, ind = '',0
  fh.close()
  Write_out(out, primer_tag_file_count)
  out, ind = '',0
  print "\tTotal sequences:\t"+str(total)+"\tNumber of REV primers found:\t"+str(pass_r)+"\tNumber of FOR primers found:\t"+str(pass_f)+"\tPassed sequences:\t"+str(pass_all)
  return()

def Get_consensus_sequence(u_seq, u_freq,tmp_file,threshold):
  out=''
  for i in range(0,len(u_seq)):
    out=out+">"+str(i)+"\n"+u_seq[i]+"\n"
  fh=open(tmp_file+"txt", "w")
  fh.write(out)
  fh.close()
  insert = ''
  if( len(u_seq) > 2000):insert = '--parttree'
  command1 = "/software/pubseq/bin/mafft-6.857/bin/mafft --retree 2 "+insert+" "+tmp_file+"txt > "+tmp_file+"aligned"
  commands.getoutput(command1)
  fh=open(tmp_file+"aligned","r")
  max_seqs={}
  pfh= Check_fasta_not_empty(fh)
  if(pfh==1):
    for header,sequence in fasta_iterator(fh):
      max_seqs[sequence.upper()]=int(header)
    fh.close()
    bases=["A","T","G","C","-"]
    base_dict = {}
    for b in range(0,len(bases)):
      base_dict[bases[b]]=b
    consensus = ''
    start = 0
    for i in range(0,len(sequence)):
      f = [0]*len(bases)
      for s in max_seqs:
        f[base_dict[s[i]]]=f[base_dict[s[i]]]+u_freq[max_seqs[s]]
      if(f[4]==0):start==1
      if(max(f)*1.0/sum(f) >= threshold):
        if(bases[f.index(max(f))]!="-"):
          consensus=consensus+bases[f.index(max(f))]
      else:
        if(start==1):
          for j in range(0,5):
            if(f[j]!=0):
              consensus=consensus+"|"+bases[j]+":"+str('%s' % float('%.3g' % (f[j]*1.0/sum(f))))
          consensus=consensus+"_"
        else:
          f = f[0:4]
          if(bases[f.index(max(f))]!="-"):
            consensus=consensus+bases[f.index(max(f))]
          else:
            for j in range(0,4):
              if(f[j]!=0):
                consensus=consensus+"|"+bases[j]+":"+str('%s' % float('%.3g' % (f[j]*1.0/sum(f))))
              consensus=consensus+"_"
    else:
      print "ERROR", u_freq, u_seq
  else:consensus =""
  return(consensus)

def Trim_sequences_BCR_TCR(Tmp_file,Fail_file,Output_trim, gene,paired,species,primer_file,primer_tag_file,tmp_file,primer_tag_file_count,sample,ref_const,reverse_primer_group,cd_hit_directory,other,barcode_group):
  (rc)= Init_rc()
  forward, reverse,  barcoded_j, barcoded_v,v_ref = Get_primers_split(rc, primer_file)
  fh_out = open(Output_trim,"w")
  fh_out.close()
  fh_out = open(Fail_file,"w")
  fh_out.close()
  threshold_BARCODE = 0.80 ### Certainty for accepting a barcode
  print barcoded_j, barcoded_v
  if(barcoded_j==1 and barcoded_v==0):
    inside = 1
    print "J barcoded"
    Single_J_barcoded_trimming_clustered(forward, reverse,  barcoded_j, barcoded_v, rc,Tmp_file,Fail_file,Output_trim,primer_tag_file,tmp_file,gene,paired,species,primer_file,primer_tag_file_count,threshold_BARCODE,ref_const,inside,v_ref,cd_hit_directory)
  if(barcoded_j==0 and barcoded_v==0):
    print "Unbarcoded"
    inside = 0
    Non_barcoded_trimming(forward, reverse,  barcoded_j, barcoded_v, rc,Tmp_file,Fail_file,Output_trim,primer_tag_file,tmp_file,gene,paired,species,primer_file,primer_tag_file_count,threshold_BARCODE,ref_const,inside,v_ref)
  if(barcoded_j==1 and barcoded_v==1):
    print "V-J barcoded"
    Double_barcoded_trimming(forward, reverse,  barcoded_j, barcoded_v, rc,Tmp_file,Fail_file,Output_trim,primer_tag_file,tmp_file,gene,paired,species,primer_file,primer_tag_file_count,threshold_BARCODE,ref_const)
  return()

def Double_barcoded_trimming(forward, reverse,  barcoded_j, barcoded_v, rc,Tmp_file,Fail_file,Output_trim,primer_tag_file,tmp_file,gene,paired,species,primer_file,primer_tag_file_count,threshold_BARCODE,ref_const):
  Read_untrimmed_file_double(rc,Tmp_file,Fail_file,Output_trim, gene,paired,species,primer_file,primer_tag_file,tmp_file,primer_tag_file_count,forward, reverse)
  Check_barcodes_MALBAC(primer_tag_file_count,primer_tag_file,Fail_file,Output_trim,threshold_BARCODE,cd_hit_directory,Tmp_file)
  inside = 0
  Print_trimmed_sequences(Output_trim, primer_tag_file, primer_tag_file_count,inside,ref_const)
  return()

def Print_trimmed_sequences(Output_trim, primer_tag_file, primer_tag_file_count,inside,ref_const):
  constant_region,regions,uniq_bcs = {},{},{}
  print "inside",inside
  if(inside in [0,1]):
    fh=open(primer_tag_file_count,"r")
    for l in fh:
      if(l[0]!="#"):
        l=l.strip().split()
        const_reg = l[len(l)-2].replace("REVERSE","REV")
        if(const_reg.count("VH")==0):
          bc = l[1]+":"+l[2]
          constant_region[l[0].split(":")[0]] = const_reg
          regions[const_reg]=1
          uniq_bcs[bc] = 1
    fh.close()
    print "Number of unique BCs:",len(uniq_bcs)
    del uniq_bcs
  fh=open(primer_tag_file,"r")
  unique_seqs,seqs = Tree(),{}
  passed,failed,uniq_failed = 0,0,{}
  if(inside in [0,1]):
    for l in fh:
      if(l[0]!="#"):
        l=l.strip().split()
        if(l[6]=="YES"):
          passed = passed+1
          consensus = l[7]
          if(consensus.count("human_BC_REV_CONST")!=0):
            print l
          seqs[l[3]+":"+l[4]] = [l[0],consensus]
        else:
          failed = failed+1
          if(l[3] in uniq_failed):uniq_failed[l[3]]= uniq_failed[l[3]]+1
          else:uniq_failed[l[3]] = 1
  print "Passed:", passed, "Failed:",failed,"(",failed*100/(failed+passed),"%)"
  #n_seqs_associated_with_failed_BCs = []
  #for n in uniq_failed:
  #  n_seqs_associated_with_failed_BCs.append(uniq_failed[n])
  #n_seqs_associated_with_failed_BCs.sort()
  #print "\t\tMean failed BC multiplicity:",mean(n_seqs_associated_with_failed_BCs),"Max failed BC multiplicity:", max(n_seqs_associated_with_failed_BCs), "Sum:",sum(n_seqs_associated_with_failed_BCs)
  fh.close()
  if(inside==0):
    for bc in seqs:
      id = seqs[bc][0].split(":")[0].split("__")[0]
      unique_seqs[seqs[bc][1]][id].value = 1
    del seqs
    out,ind = '',0
    for seq in unique_seqs:
      f = len(unique_seqs[seq])
      for id in unique_seqs[seq]:
        break
      out=out+">"+id.split("__")[0]+"__"+str(f)+"\n"+seq+"\n"
      ind = ind+1
      if(ind>100):
        Write_out(out, Output_trim)
        out, ind = '',0
    Write_out(out, Output_trim)
  elif(inside==1):
    print "NOT YET DEVELOPED!"
    fh=open(ref_const,"r")
    word_size,ref_words = 10,{}
    for header,seq in fasta_iterator(fh):
      words,seq = [],seq.upper()
      for i in range(word_size,len(seq)):
        words.append(seq[i-word_size:i])
      ref_words[header.split("|")[0]]=words
    for bc in seqs:
      seq,id = seqs[bc][1], seqs[bc][0]
      score = []
      #print bc, seqs[bc]
      ##id = seqs[bc][0].split(":")[0].split("__")[0]
      #if(id in constant_region):
      #  unique_seqs[seqs[bc][1]][constant_region[id]][id].value = 1
        #print seqs[bc][1], constant_region[id], id
    del seqs

def Read_untrimmed_file_double(rc,Tmp_file,Fail_file,Output_trim, gene,paired,species,primer_file,primer_tag_file,tmp_file,primer_tag_file_count,forward, reverse):
  for f in [primer_tag_file_count]:
    fh=open(f, "w")
    fh.close()
  fh = open (Tmp_file, "r")
  seqs = Tree() 
  maxl=1000000
  if(gene=="HEAVY" or gene == "IGH"):minl= 120
  if(gene=="KAPPA" or gene=="IGK"):minl= 110 
  if(gene=="LAMBDA" or gene == "IGL"):(minl, maxl)= (90, 150)
  J_found,v_found,tot,indexing=0,0,0,0
  seqs1,t=Tree(),0
  for header,seq in fasta_iterator(fh):
    seq=seq.upper()
    seqs1[seq][header].value=1
    t=t+1
  out,ind = "#ID\tJ_tag\tV_tag\tSequence\ttype_rev\ttype_for\tprimer_rev\tprimer_for\n", 0
  total,pass_r, pass_f, pass_all = 0,0,0,0
  for seq in seqs1:
    for header in seqs1[seq]:
      break
    header = header+":"+str(len(seqs1[seq]))
    number = len(seqs1[seq])
    passes,j_tag,v_tag=0,'',''
    type_rev, type_for = '',''
    primer_rev,primer_for = '',''
    total = total+1
    for i in range(0,len(reverse)):
      pj=Get_match(reverse[i][2], seq)
      if(max(pj+[-1])==-1 or len(pj)==0):
        seq = Reverse_comp(seq, rc) 
        pj=Get_match(reverse[i][2], seq)
      if(max(pj+[-1])!=-1):
        bc_len = len(reverse[i][3])+3#len(reverse[i][2])/5
        if(pj[0]>bc_len):
          #start = pj[0]
          #primer_match_seq = seq[start:start+len(reverse[i][2])]
          j_tag = seq[pj[0]+len(reverse[i][2]):pj[0]+len(reverse[i][2])+bc_len]
          #print primer_match_seq, reverse[i][2], seq[start:start+len(reverse[i][2])+30],j_tag
          seq = seq[0:pj[0]+bc_len]
          type_rev,primer_rev = reverse[i][1], reverse[i][5]
          passes = 1
          break
    #if(passes!=1):
    if(passes ==1000):
      for i in range(0,len(reverse)):
        words =reverse[i][7]
        p=[]
        for w in words:
          pj=Get_match(w[0], seq)
          if(max(pj+[-1])!=-1 and len(pj)>0):
            if(pj[0]<len(seq)/2):p=p+[max([0,pj[0]-w[1]])]
            else:pj[0]=-1
        if(len(p)>1):
          bc_len = len(reverse[i][3])+3
          pj = max(set(p), key=p.count)
          start = pj
          j_tag = seq[pj-bc_len:pj+1]
          #start = pj-bc_len
          primer_match_seq = seq[start:start+len(reverse[i][6])]
          print primer_match_seq, reverse[i][2], seq[start:start+len(reverse[i][6])+30]
          seq = seq[pj+bc_len:len(seq)]
          seq = Reverse_comp(seq, rc)
          type_rev,primer_rev = reverse[i][1], reverse[i][5]
          passes = 1
          break
      if(passes!=1):
        seq = Reverse_comp(seq, rc) 
        for i in range(0,len(reverse)):
          words =reverse[i][6]
          p=[]
          for w in words:
            pj=Get_match(w[0], seq)
            if(max(pj+[-1])!=-1 and len(pj)>0):
              if(pj[0]<len(seq)/2):p=p+[max([0,pj[0]-w[1]])]
              else:pj[0]=-1
          if(len(p)>1):
            bc_len = len(reverse[i][3])+3
            pj = max(set(p), key=p.count)
            j_tag = seq[pj-bc_len:pj+1]
            start = pj-bc_len 
            primer_match_seq = seq[start:start+len(reverse[i][6])]
            print primer_match_seq, reverse[i][2], seq[start:start+len(reverse[i][6])+30]
            seq = seq[pj+bc_len:len(seq)]
            seq = Reverse_comp(seq, rc)
            type_rev,primer_rev = reverse[i][1], reverse[i][5]
            passes = 1
            break
    if(passes==1):
      pass_r = pass_r+1
      for i in range(0,len(forward)):
        pv=Get_match(forward[i][0], seq)
        if(max(pv+[-1])!=-1):
          v_tag='-'
          seq = seq[pv[0]+len(forward[i][0])-1 :len(seq)]
          type_for, primer_for = forward[i][1],forward[i][2]
          passes = 2
          break
      if(passes!=2):
        for i in range(0,len(forward)):
          words =forward[i][6]
          p=[] 
          for w in words:
            pv=Get_match(w[0], seq) 
            if(max(pv+[-1])!=-1 and len(pv)>0):
              if(pv[0]<len(seq)/2):p=p+[max([0,pv[0]-w[1]])]
              else:pv[0]=-1
          if(len(p)>3):
            pv = max(set(p), key=p.count)
            bc_len = len(forward[i][3])+1
            v_tag = seq[pv-bc_len:pv]
            type_for, primer_for = forward[i][1],forward[i][2]
            #print type_for, primer_for,v_tag,seq[pv-bc_len-6:pv+4], p
            seq = seq[pv+len(forward[i][2])-1 :len(seq)]
            type_for, primer_for = forward[i][1],forward[i][2]
            passes = 2
            break
    if (passes ==2):
      pass_f = pass_f+1
      if(len(seq)>minl and len(seq)<maxl):
        if(seq.count("N")==0):
          #if(type_rev==type_for):
          out=out+"\t".join(map(str,[header, j_tag, v_tag, seq, type_rev,type_for,primer_rev,primer_for]))+"\n"
          pass_all = pass_all+1
          ind = ind+1
          if(ind>200):
            Write_out(out, primer_tag_file_count)
            out, ind = '',0
  fh.close()
  Write_out(out, primer_tag_file_count)
  out, ind = '',0
  print "\tTotal sequences:\t"+str(total)+"\tNumber of REV primers found:\t"+str(pass_r)+"\tNumber of FOR primers found:\t"+str(pass_f)+"\tPassed sequences:\t"+str(pass_all)
  return()



def Non_barcoded_trimming(forward, reverse,  barcoded_j, barcoded_v, rc,Tmp_file,Fail_file,Output_trim,primer_tag_file,tmp_file,gene,paired,species,primer_file,primer_tag_file_count,threshold_BARCODE,ref_const,inside,v_ref):
  for f in [primer_tag_file_count, Output_trim]:
    fh=open(f, "w")
    fh.close()
  fh = open (Tmp_file, "r")
  seqs = Tree() 
  minl,maxl=110,1000000
  if(gene=="HEAVY" or gene == "IGH"):minl= 120
  if(gene=="KAPPA" or gene=="IGK"):minl= 110 
  if(gene=="LAMBDA" or gene == "IGL"):(minl, maxl)= (90, 150)
  J_found,v_found,tot,indexing=0,0,0,0
  seqs1,t=Tree(),0
  for header,seq in fasta_iterator(fh):
    seq=seq.upper()
    seqs1[seq][header].value=1
    t=t+1
  out,ind = "#ID\tJ_tag\tV_tag\tSequence\ttype_rev\ttype_for\tprimer_rev\tprimer_for\n", 0
  total,pass_r, pass_f, pass_all = 0,0,0,0
  seqs_new = Tree()
  for seq in seqs1:
    for header in seqs1[seq]:
      break
    header = header+":"+str(len(seqs1[seq]))
    number = len(seqs1[seq])
    passes,j_tag,v_tag=0,'',''
    type_rev, type_for = '',''
    primer_rev,primer_for = '',''
    total = total+1
    for i in range(0,len(reverse)):
      pj=Get_match(reverse[i][2], seq)
      if(max(pj+[-1])==-1 or len(pj)==0):
        seq = Reverse_comp(seq, rc) 
        pj=Get_match(reverse[i][2], seq)
      if(max(pj+[-1])!=-1):
        if(pj[0]>0):
          j_tag = '-'
          #seq = seq[pj[0]:len(seq)]
          seq = seq[0:pj[0]]
          type_rev,primer_rev = reverse[i][1], reverse[i][3]
          passes = 1
          break
    if(passes!=1):
      for i in range(0,len(reverse)):
        words =reverse[i][4]
        p=[]
        for w in words:
          pj=Get_match(w[0], seq)
          if(max(pj+[-1])!=-1 and len(pj)>0):
            if(pj[0]<len(seq)/2):p=p+[max([0,pj[0]-w[1]])]
            else:pj[0]=-1
        if(len(p)>1):
          pj = max(set(p), key=p.count)
          if(pj>0):
            j_tag = '-'
            seq = seq[pj:len(seq)]
            seq = Reverse_comp(seq, rc)
            type_rev,primer_rev = reverse[i][1], reverse[i][3]
            passes = 1
            break
      if(passes!=1):
        seq = Reverse_comp(seq, rc) 
        for i in range(0,len(reverse)):
          words =reverse[i][4]
          p=[]
          for w in words:
            pj=Get_match(w[0], seq)
            if(max(pj+[-1])!=-1 and len(pj)>0):
              if(pj[0]<len(seq)/2):p=p+[max([0,pj[0]-w[1]])]
              else:pj[0]=-1
          if(len(p)>1):
            pj = max(set(p), key=p.count)
            if(pj>0):
              j_tag = '-'
              seq = seq[pj:len(seq)]
              seq = Reverse_comp(seq, rc)
              type_rev,primer_rev = reverse[i][1], reverse[i][3]
              passes = 1 
              break
    if(passes==1):
      pass_r = pass_r+1
      for i in range(0,len(forward)):
        pv=Get_match(forward[i][0], seq)
        if(max(pv+[-1])!=-1):
          v_tag='-'
          seq = seq[pv[0]+len(forward[i][0])-1 :len(seq)]
          type_for, primer_for = forward[i][1],forward[i][2]
          passes = 2
          break
      if(passes!=2):
        for i in range(0,len(forward)):
          words =forward[i][3]
          p=0 
          for w in words:
            pv=Get_match(w[0], seq) 
            if(max(pv+[-1])!=-1 and len(pv)>0):
              if(pv[0]<len(seq)/2):p=p+1
              else:pv[0]=-1
          if(max(pv+[-1])!=-1 and len(pv)>0 and p>1):
            v_tag='-'
            seq = seq[pv[0]+len(forward[i][2])-1 :len(seq)]
            type_for, primer_for = forward[i][1],forward[i][2]
            passes = 2
            break
      if (passes!=2):
        p=[]
        for i in range(0,len(v_ref)):
          words= v_ref[i][3]
          p1 = []
          for w in words:
            pv=Get_match(w[0], seq)
            if(max(pv+[-1])!=-1 and len(pv)>0):
              if(pv[0]<len(seq)/2):p=p+[pv[0]]
          if(len(p1)>0):p1 = p1+[min(p1)]
        if(len(p)>=2):
          pv = min(p)
          v_tag='-'
          seq = seq[pv:len(seq)]
          type_for, primer_for = v_ref[i][1],v_ref[i][2]
          passes = 2
    if (passes ==2):
      pass_f = pass_f+1
      #print seq
      if(len(seq)>minl and len(seq)<maxl):
        if(seq.count("N")==0):
          out=out+"\t".join(map(str,[header, j_tag, v_tag, seq, type_rev,type_for,primer_rev,primer_for]))+"\n"
          seqs_new[seq][header].value = 1
          pass_all = pass_all+1
          ind = ind+1
          if(ind>200):
            Write_out(out, primer_tag_file_count)
            out, ind = '',0
  fh.close()
  Write_out(out, primer_tag_file_count)
  out, ind = '',0
  print "\tTotal sequences:\t"+str(total)+"\tNumber of REV primers found:\t"+str(pass_r)+"\tNumber of FOR primers found:\t"+str(pass_f)+"\tPassed sequences:\t"+str(pass_all)
  fh=open(primer_tag_file_count,"r")
  out, ind = '',0
  for seq in seqs_new:
    f = 0
    for id in seqs_new[seq]:
      f = f+int(id.split(":")[1])
    out = out+">"+id.split(":")[0]+"__"+str(f)+"\n"+seq+"\n"
    ind = ind+1
    if(ind>200):
      Write_out(out, Output_trim)
      out, ind = '',0
  Write_out(out, Output_trim)
  out, ind = '',0
  return()

def Get_primers_split(rc, primer_file):
  fh=open(primer_file, "r")
  forward, reverse,v_ref = [],[],[]
  barcoded_j, barcoded_v = 0,0
  word_size = 8
  for l in fh:
    if(l[0]!="#"):
      l=l.strip().split()
      if(len(l)>2):
        header, sequence, uni, BC,gene_specific =l[0],l[1],l[2],l[3],l[4]
        gene_specific=gene_specific.upper()
        words,words1 = [],[]
        for i in range(word_size,len(gene_specific)):
          words.append([gene_specific[i-word_size:i],i-word_size])
        for i in range(word_size,len(uni)):
          words1.append([uni[i-word_size:i],i-word_size])
        if(header.count("J")!=0 or header.count("REV_CONST")!=0 or header.count("REVERSE")!=0):
          if(header.count("REV_CONST")!=0):inside=1
          if(sequence.count("N")!=0):barcoded_j=1
          sequence=sequence.upper()
          if(header.count("HOUSEKEEPING")!=0):clas = "HOUSEKEEPING"
          else:clas = "IMMUME_REC"
          l=len(gene_specific)
          reverse =reverse+[[sequence, clas, uni, BC,gene_specific,header,words,words1]]
        else:
          if(sequence.count("N")!=0):barcoded_v = 1
          sequence=sequence.upper()
          if(header.count("HOUSEKEEPING")!=0):clas = "HOUSEKEEPING"
          else:clas = "IMMUNE_REC"
          l=len(gene_specific)
          forward = forward+[[sequence, clas, uni, BC,gene_specific,header,words,words1]]
      else:
        header, sequence=l[0],l[1]
        sequence=sequence.upper()
        words = []
        for i in range(word_size,len(sequence)):
          words.append([sequence[i-word_size:i],i-word_size])
        if(header.count("J")!=0 or header.count("REV_CONST")!=0 or header.count("REVERSE")!=0):
          if(header.count("CONST")!=0):inside=1
          if(sequence.count("N")!=0):barcoded_j=1
          sequence=sequence.upper()
          clas = "IMMUME_REC"
          reverse =reverse+[[sequence,clas, sequence,header, words]]
        elif(header.count("REF")==0):
          if(sequence.count("N")!=0):barcoded_v = 1
          sequence=sequence.upper()
          clas = "IMMUNE_REC"
          forward = forward+[[sequence, clas,header,words]]
        else:
          sequence=sequence.upper()
          clas = "IMMUNE_REC"
          v_ref = v_ref+[[sequence, clas,header,words]]
  fh.close()
  return(forward, reverse,  barcoded_j, barcoded_v,v_ref)

def Check_location_of_primer_binding(primer_file):
  fh=open(primer_file,"r")
  mode = "WITHIN"
  for header,sequence in fasta_iterator(fh):
    if(header.count("REV_CONST")!=0):
      mode = "OUTSIDE"
  fh.close()
  return(mode)

def Blast_match_J_const(out, seqs, Trim1, Trim2, refj,e_value,indent):
  fh=open(Trim1+"_blast_J","w")
  fh.write(out)
  fh.close() 
  command1 = "blastall -p blastn -a 10 -d "+refj+" -e "+str(e_value)+" -i "+Trim1+"_blast_J -o "+Trim1+"_blast_J_results -b 1 -m 8"
  commands.getoutput(command1)
  fh = open(Trim1+"_blast_J_results","r")
  out= ''
  for l in fh:
    l=l.strip().split()
    header, gene = l[0].split("|")[1].split("_"), map(int, l[0].split("|")[0].split("__")[1].split("_"))
    gene = header[gene.index(max(gene))]
    if(l[1].count(gene[0:3])!=0):
      if(int(l[6])<int(l[7]) and int(l[8])<int(l[9])):
        j_end = int(l[7])
        seq_j = seqs[l[0]][0]
        seq_j = seq_j[0:j_end+seqs[l[0]][1]]
        out = out+">"+l[0]+"\n"+seq_j+"\n"
  fh.close()
  Write_out(out, Trim2)
  return()

def Filter_IgJ_genes(Trim1, Trim2, refj,primer_file,ref_const,primer_tag_file_count,blastall,refv,primer_length):
  mode ="WITHIN" 
  fh=open(Trim2,"w")
  fh.close()
  if(1==1):
    if(mode == "WITHIN"):
      out,ind,batch,batch_size = '',0,0,500
      seqs,indent={},120
      fh=open(Trim1, "r")
      e_value = 10
      c = 0
      for header,sequence in fasta_iterator(fh):
        #inf = len(sequence)-indent
        inf = 0
        #if(inf<0):10
        out=out+">"+header+"\n"+sequence[inf:len(sequence)].replace("...","")+"\n"
        seqs[header]=sequence
        ind,batch = ind+1,batch+1
        c = c+1
        if(batch>=batch_size):
          Blast_match_J(out, seqs, Trim1, Trim2, refj,e_value,blastall,refv,primer_length)
          out,batch,seqs = '',0,{}
          print "BATCH"
      fh.close()
      if(len(out)>2):
        Blast_match_J(out, seqs, Trim1, Trim2, refj,e_value,blastall,refv,primer_length)
        out,batch = '',0
        print "BATCHEND"
    elif(mode=="OUTSIDE"):
      out,ind,batch,batch_size = '',0,0,500
      seqs,indent={},150
      fh=open(Trim1, "r")
      e_value = 0.1
      for header,sequence in fasta_iterator(fh):
        inf = len(sequence)-indent
        if(inf<0):10
        out=out+">"+header+"\n"+sequence[inf:len(sequence)]+"\n"
        seqs[header]=[sequence,inf]
        ind,batch = ind+1,batch+1
        if(batch>=batch_size):
          Blast_match_J_const(out, seqs, Trim1, Trim2, refj,e_value,indent)
          out,batch,seqs = '',0,{}
      fh.close()
      if(len(seqs)>0):
        Blast_match_J_const(out, seqs, Trim1, Trim2, refj,e_value,indent)
        out,batch = '',0
  return()

def Check_fasta_not_empty(fh):
  pfh = 0
  for l in fh:
    if(len(l)!=0):pfh =1
    break
  return(pfh)

def Blast_match_J(out, seqs, Trim1, Trim2, refj,e_value, blastall,refv,primer_length):
  fh=open(Trim1+"_blast_J","w")
  fh.write(out)
  fh.close()
  command1 =  blastall+ " -db "+refj+" -evalue "+str(e_value)+" -query "+Trim1+"_blast_J -out  "+Trim1+"_blast_J_results -max_target_seqs 1 -outfmt 6 -word_size 4"
  commands.getoutput(command1)
  command2 =  blastall+ " -db "+refv+" -evalue "+str(e_value)+" -query "+Trim1+"_blast_J -out  "+Trim1+"_blast_V_results -max_target_seqs 1 -outfmt 6 -word_size 8"
  commands.getoutput(command2)
  fh = open(Trim1+"_blast_V_results","r")
  done={}
  for l in fh:
    l=l.strip().split()
    if(l[0] not in done):
      print int(l[3])
      if(int(l[3])>50):
        start = int(l[5])+primer_length
        done[l[0]]= start
  fh.close()
  fh = open(Trim1+"_blast_J_results","r")
  out,done1='',{}
  for l in fh:
    l=l.strip().split()
    if(l[0] in done and l[0] not in done1):
      start = done[l[0]]
      seq = seqs[l[0]]
      out=out+">"+l[0]+"\n"+seq[start:len(seq)]+"\n"
      done1[l[0]]= 1
  fh.close()
  fh=open(Trim2,"a")
  fh.write(out)
  fh.close()
  return()

def Blast_match(V_region, ref, tmp_file,v_match):
  fh=open(tmp_file, "w")
  fh.write(V_region)
  fh.close()
  command1 = "blastall -p blastn -a 10 -d "+ref+" -e 1 -i "+tmp_file+" -o "+tmp_file+"_blast -b 1 -m 8"
  commands.getoutput(command1)
  fh = open(tmp_file+"_blast","r")
  passes={}
  for l in fh:
    l=l.strip().split()
    if(int(l[3])>v_match):
      passes[l[0]]= l[1]
  return(passes)

def Check_type_of_priming(primer_file):
  fh=open(primer_file,"r")
  CONST,BARCODE = "FALSE","FALSE"
  for l in fh:
    l=l.strip().split()
    if(l[0].count("_CONST")!=0):
      CONST = "TRUE"
    if(len(l)>2):
      BARCODE="TRUE"
  fh.close()
  return(CONST,BARCODE)

def Reduce_sequences(Trim2, Trim3,primer_file):
  CONST,BARCODE = Check_type_of_priming(primer_file)
  min_length = 300
  min_length = 180
  if(CONST=="TRUE"):
    lengths = {}
    fh=open(Trim3,"w")
    fh.close()
    fh=open(Trim2,"r")
    seqs,head = Tree(),''
    for header,seq in fasta_iterator(fh):
      if(header.count(">")!=0):
        header = header.split(">")
        header = header[len(header)-1]
        print header
      seq=seq.upper()
      l = len(seq)
      if(l in lengths):lengths[l] = lengths[l]+1
      else:lengths[l] = 1
      if(len(seq)>=min_length):
        seqs[seq][header].value = 1
        head = header.split("|")[1]
    fh.close()
    out,ind,times='',0,len(head.split("_"))
    print "number of chains", times 
    for seq in seqs:
      f = [0]*times
      for id in seqs[seq]:
        f1 = map(int, id.split("__")[1].split("|")[0].split("_"))
        f = [a + b for a, b in zip(f,f1)]
        #f = map(add, f, f1)
      header = ">"+id.split("__")[0]+"__"+"_".join(map(str,f))+"|"+head
      out=out+header+"\n"+seq+"\n"
      ind = ind+1
      if(ind>500):
        Write_out(out, Trim3)
        out, ind = '',0
    Write_out(out, Trim3)
    for l in lengths:
      print l,"\t",lengths[l]
  else:
    command1 = "cp "+Trim2+" "+Trim3
    commands.getoutput(command1)
  return()

def Get_sequences_ref(file,word):
  fh = open(file,"r")
  dict=Tree()
  for header,sequence in fasta_iterator(fh):
    for i in range(1, len(sequence)-word-1):
      dict[sequence[i:(i+word)]][header].value = 1
  fh.close()
  return(dict)

def Get_max_ORF(ORF,word,dict):
  score=[] 
  codon=[]
  if(len(ORF)>1):
    for si in ORF:
      s=si[0]
      sc=0
      for i in range(1, len(s)-word):
        if(s[i:(i+word)] in dict):sc = sc+1
      score.append(sc)
      codon.append(si[1])
    maxi=score.index(max(score))
    orf=ORF[maxi][0]
    mscore=max(score)
    codon=ORF[maxi][1]
  else:
    for si in ORF:
      orf=si[0]
      mscore=1000
      codon=si[1]
  return(orf, mscore, codon)

def Translate(seq,codon):
  p_seq=""
  for cod in range(0, len(seq)/3-1):
    cod = cod*3
    c1 = seq[cod+0]+seq[cod+1]+seq[cod+2]
    if(c1 in codon):
      p_seq=p_seq+str(codon[c1])
  return(p_seq)

def Calculate_ORF_Length (codon,sequence,type,gene,read):
  three_frames = [ sequence[i:] for i in range(3) ]
  three_frames_translated = [Translate(frame,codon) for frame in three_frames]
  accepted=[]
  found=0
  for i in range(0,len(three_frames_translated)):
    frame=three_frames_translated[i]
    if(frame.count("-")==0):
      accepted.append((frame,i))
  accepts = 1
  if(len(accepted)==0):
    accepts = 0
  return(accepted,accepts)

def Get_score(a1,a2):
  score = 0
  for i in range(0,len(a1)):
    if(a1[i]!="-"):
      if(a1[i]==a2[i]):
        score=score+1
  return(score)

def Get_freq (id):
  id = id.split("__")
  return(int(id[1]))

def Get_nucleotide_sequences(Output_trim, Filtered_out1,nn_orf_filtered,dir_ind,gene,ref,refj, ref_protein,refjp,tmp_file):
  fh=open(nn_orf_filtered,"w")
  fh.close()
  fh=open(Filtered_out1,"r")
  ids = {}
  for header,sequence in fasta_iterator(fh):
    ids[header] = 1
  fh.close()
  out, ind = '',0
  fh=open(Output_trim,"r")
  done = {}
  for header,sequence in fasta_iterator(fh):
    if(header in ids):
      out=out+">"+header+"\n"+sequence+"\n"
      done[header] = 1
      ind = ind+1
      if(ind>500):
        Write_out(out, nn_orf_filtered)
        out, ind = '',0
  Write_out(out, nn_orf_filtered)
  fh.close()
  print len(ids), len(done)
  return()

def Get_codons (codon_table):
  fh = open(codon_table,"r")
  codon = {}
  for l in fh:
    l=l.strip()
    l1 = l.split()
    l2 = list(l1[0])
    codon[l2[0]+l2[1]+l2[2]]=l1[1]
  fh.close()
  return (codon)

def Count_isotypes(file,sample):
  fh = open(file)
  chains = {}
  total = 0
  for header,seq in fasta_iterator(fh):
    freq = map(int, header.split("__")[1].split("|")[0].split("_"))
    iso = header.split("|")[1].split("_")
    nz = [i for i in range(len(freq)) if freq[i]!=0]
    total = total+sum(freq)
    for i in nz:
      if(iso[i] in chains):chains[iso[i]]=chains[iso[i]]+freq[i]
      else:chains[iso[i]]=freq[i]
  t = 0
  for c in chains: 
    print sample,"\t",c,"\t",chains[c],"\t",chains[c]*100.0/total
    t = t+chains[c]
  print t, total
  return()

def Get_protein_sequences(Output_trim, Filtered_out1,nn_orf_filtered,dir_ind,gene,ref,refj, ref_protein,refjp,tmp_file,codon_table):
  fh = open(Filtered_out1,"w")
  fh.close()
  fh=open(nn_orf_filtered, "w")
  fh.close()
  (codon)=Get_codons (codon_table)
  (word,v_match) = (4,45)
  (dict1)=Get_sequences_ref(ref_protein,word)
  fh = open(Output_trim,"r")
  number_seqs,number_accepted, batch, indb, V_region, seqs, ORFa, ORFb, batch_number, orf_fail, p_seq=0,0,1000,0,'',Tree(),{},{},0,0,{}
  out = ''
  chains = {}
  for header,seq in fasta_iterator(fh):
    if(seq.count("N")==0):
      seq=seq.replace("-","")
      header= header.split("#")[0]
      if(header.count("__")==0):
        freq = Get_freq(header)
        header = header.replace(":","")+"__"+str(freq)
      (ORF1,accept1)=Calculate_ORF_Length (codon, seq, "V",gene,1)
      ORFa[header] = ORF1
      number_seqs=number_seqs+1#freq
      if(accept1==0):
        orf_fail = orf_fail+1
      if(accept1 !=0 ):
        if(len(ORF1)>1):
          min_score,found = 2,0
          (O1,score1,codon1) =  Get_max_ORF(ORF1,word,dict1)
          if(score1>min_score):
            min_score = score1
            p_seq[header]=O1
            seq = seq[codon1:len(seq)]
            found = 1
          if(found==1):number_accepted = number_accepted+1#int(header.split("__")[1])
        else:
          p_seq[header]=ORF1[0][0]
          number_accepted = number_accepted+1#int(header.split("__")[1])
          seq = seq[ORF1[0][1]:len(seq)]
        seqs[seq][header].value=1
        indb=indb+1
        if(indb >= batch):
          batch_number = batch_number+1
          out=''
          for id in p_seq:
            out=out+">"+id+"\n"+p_seq[id]+"\n"
          Write_out(out, Filtered_out1)
          p_seq={}
          out=''
  fh.close()
  Write_out(out, Filtered_out1)
  out,ind='',0
  for id in p_seq:
    out=out+">"+id+"\n"+p_seq[id]+"\n"
    ind = ind+1
    if(ind>500):
      Write_out(out, Filtered_out1)
      out, ind = '',0
  Write_out(out, Filtered_out1)
  return()

def ORF_calculation_single (Output_trim, Filtered_out1,nn_orf_filtered,dir_ind,gene,ref,refj, ref_protein,refjp,tmp_file,codon_table,ORF_filter):
  if(ORF_filter =="True"):
    Get_protein_sequences(Output_trim, Filtered_out1,nn_orf_filtered,dir_ind,gene,ref,refj, ref_protein,refjp,tmp_file,codon_table)
    Get_nucleotide_sequences(Output_trim, Filtered_out1,nn_orf_filtered,dir_ind,gene,ref,refj, ref_protein,refjp,tmp_file)
  else:
    command = ("cp "+Output_trim+" "+nn_orf_filtered)
    print command
    commands.getoutput(command)
  return()

def Get_number_sequences_unique(file):
  c = commands.getoutput("wc -l "+file).split()
  if(c[0] not in ['wc:',"ls:"]):
    fh=open(file,"r")
    seqs = {}
    for header,sequence in fasta_iterator(fh):
      seqs[header] = 1
    fh.close()
    count = len(seqs)
  else:
    count = -1
  return(count)

def Get_number_sequences(file):
  count = 0
  c = commands.getoutput("wc -l "+file).split()
  if(c[0] not in ['wc:',"ls:"]):
    fh=open(file,"r")
    for header,sequence in fasta_iterator(fh):
      count = count+1
    fh.close()
  else:count = -1
  return(count)

def Get_reduced_number_sequences_multi_constants(file):
  count = 0
  c = commands.getoutput("wc -l "+file).split()
  if(c[0]!="ls:" and c[0]!="0"):
    fh=open(file,"r")
    for header,sequence in fasta_iterator(fh):
      count=count+ sum(map(int,header.split("__")[1].split("|")[0].split("_")))
    fh.close()
  return(count)

def Count_barcodes(primer_tag_file):
  print primer_tag_file
  fh=open(primer_tag_file, 'r')
  bc_uniq, bc_count,bc_failed = {},0,0
  uniq_seq = 0
  for l in fh:
    if(l[0]!="#"):
      l=l.strip().split()
      bc = l[3]+":"+l[4]
      if(l[1]!="NA"):
        freq = int(l[2])
        uniq_seq = uniq_seq+freq
        bc_uniq[bc] = 1
        bc_count = bc_count+freq
        if(l[6]!="YES"):
          bc_failed = bc_failed+1
  fh.close()
  n_barcodes,uniq_sequences, total_sequences_included_before_bc = len(bc_uniq), uniq_seq, bc_count
  return(n_barcodes,uniq_sequences, total_sequences_included_before_bc)

def Get_read_report(Seq_file1,Seq_file2,Tmp_file,Trim1, nn_orf_filtered,filtering_report,id, species, gene, dir,primer_tag_file_count,primer_file,method,barcode_group,Trim3,other):
  if(method =="Multiplex_FORWARD_BARCODE_GROUPED"):Seq_file1,Seq_file2 = dir+"FASTQ_FILES/Sequences_"+barcode_group+"_1.fasta", dir+"FASTQ_FILES/Sequences_"+barcode_group+"_2.fasta"
  rc = {}
  CONST,BARCODE = Check_type_of_priming(primer_file)
  if(BARCODE=="TRUE"):barcoded_j,barcoded_v=1,0
  elif(BARCODE=="FALSE"):barcoded_j,barcoded_v=0,0
  if(barcoded_j+barcoded_v==0):
    raw1,joined,gene_matching,nORF=Get_number_sequences(Seq_file1),Get_number_sequences(Tmp_file), Get_reduced_number_sequences_multi_constants(Trim1), Get_reduced_number_sequences_multi_constants(Trim3)
    orf = Get_reduced_number_sequences_multi_constants(nn_orf_filtered)
    number_unique_seqs = Get_number_sequences(nn_orf_filtered)
    out = "Directory\tSample\tSpecies\tGene\t% reads retained\tN raw reads (batch)\tN joined reads\tN reads post-V marching\tN reads post-J matching\tN BCR filtered (post ORF filtering)\tN unique BCRs\n"
    out=out+dir+"\t"+id+"\t"+species+"\t"+gene+"\t"+str(orf*100.0/joined)+"\t"+str(raw1)+"\t"+str(joined)+"\t"+str(gene_matching)+"\t"+str(nORF)+"\t"+str(orf)+"\t"+str(number_unique_seqs)+"\n"
  else:
    Tmp_file1 = Tmp_file
    n_barcodes,uniq_sequences, total_sequences_included_before_bc =  Count_barcodes(primer_tag_file)
    raw1=Get_number_sequences(Seq_file1)
    joined = Get_number_sequences_unique(Tmp_file1)
    gene_matching=Get_reduced_number_sequences_multi_constants(Trim1)
    count_bc_found,count_uniq_bcs = -1,-1
    nORF=Get_reduced_number_sequences_multi_constants(Trim3)
    orf = Get_reduced_number_sequences_multi_constants(nn_orf_filtered)
    number_unique_seqs = Get_number_sequences(nn_orf_filtered)
    out = "Directory\tSample\tSpecies\tGene\t% reads retained\tN raw reads joined (batch)\tN joined reads\tN reads with UMIs\tN uniq UMIs\tN reads post-V marching\tN reads post-J matching\tN BCR filtered (post ORF filtering)\tN unique BCRs\n" 
    orf_perc = str(orf*100.0/joined)
    out=out+dir+"\t"+id+"\t"+species+"\t"+gene+"\t"+orf_perc+"\t"+str(raw1)+"\t"+str(joined)+"\t"+str(uniq_sequences)+"\t"+str(n_barcodes)+"\t"+str(gene_matching)+"\t"+str(nORF)+"\t"+str(orf)+"\t"+str(number_unique_seqs)+"\n"
  fh=open(filtering_report,"w")
  fh.write(out)
  fh.close()
  print out
  return()

def Cluster_i (Reduced_file,tmp_file,diff,cd_hit_directory):
  command= cd_hit_directory+"cd-hit -i "+Reduced_file+" -o "+tmp_file+" -c "+str(diff)+" -d 180 -T 10  -M 0 -AL 40 "
  print command
  os.system(command)
  return()

def Get_seqs_single (file):
  fh=open(file,"r")
  seqs={}
  for header,sequence in fasta_iterator(fh):
    if(header.count("M02446169000000000CK5M02446177000000000CN2DJ121141150520324")!=0):
      print header
    seqs[header.split("__")[0]]=sequence
  fh.close()
  return(seqs)

def Check_same(ids, seqs_sub, same,inv):
  l=len(ids)
  ind = 0
  for i in range(0,l):
    for j in range(i,l):
      if(seqs_sub[i]==seqs_sub[j]):
        same [ids[i]][ids[j]].value = 1
        inv[ids[j]]=ids[i]
        ind = ind+1
      else:
        s1 = seqs_sub[i]
        s2 = seqs_sub[j]
        (sa,sb,p)=Trim_sequences(s1,s2,len(s1),len(s2))
        if(sa==sb):
          same [ids[i]][ids[j]].value = 1
          inv[ids[j]]=ids[i]
          ind = ind+1
    if(ind>l/2):break
  return(same,inv)

def Trim_sequences(s1,s2,l1,l2):
  p=0
  i=15
  sample1=s1[i:i+25]
  index=s2.find(sample1)
  if (index!=-1):
    if(index>i):
      if(index-i <=20):s2=s2[index-i:l2]
    else:
      if(i-index <=20):s1=s1[i-index:l1]
    min_len=min([len(s1),len(s2)])
    if ((max([len(s1),len(s2)]) - min_len) <25):
      s1=s1[0:min_len]
      s2=s2[0:min_len]
      p=1
  else:
    i=l1-50
    sample1=s1[i:i+25]
    index=s2.find(sample1)
    if (index!=-1):
      if(index>i):
        if(index-i <=20):s2=s2[index-i:l2]
      else:
        if(i-index <=20):s1=s1[i-index:l1]
      min_len=min([len(s1),len(s2)])
      if ((max([len(s1),len(s2)]) - min_len) <25):
        s1=s1[0:min_len]
        s2=s2[0:min_len]
        p=1
      else:
        p=0
    else:
      p=0
  return (s1, s2, p)

def Get_clusters (file):
  fh=open(file,"r")
  cluster=Tree()
  for l in fh:
    l=l.strip()
    l=l.split()
    id = l[2].replace("...","")
    if(id.count(">")==2):id = id.split(">")[2]
    id = id.replace(">","")
    cluster[l[0]][id].value=1
  fh.close()
  return(cluster)

def Get_inv_clusters(file):
  fh=open(file,"r")
  inv={}
  clust=Tree()
  for l in fh:
    l=l.strip()
    l=l.split()
    id = l[2].replace("...","")
    id = id.replace(">","")
    inv[id] = l[0]
    clust[l[0]][id].value=1
  fh.close()
  return(inv)

def Merge_clusters(clust1, inv2):
  cluster=Tree()
  for c in clust1:
    for id in clust1[c]:
      if (id in inv2):
        cluster[c+"|"+inv2[id]][id].value = 1
  return(cluster)

def Get_sequence_same_single(seqs, tmp_file1,reduced_sequences, tmp_pre, read_number_division):
  fh=open(reduced_sequences, "w")
  fh.close()
  (clust)= Get_clusters(tmp_file1+".bak.clstr")
  (out, ind)=('',0)
  same = Tree()
  inv = {}
  for c in clust:
    if(len(clust[c])>1):
      ids=[]
      seqs_sub=[]
      for id in clust[c]:
        ids.append(id)
        seqs_sub.append(seqs[id.split("__")[0]])
      (same,inv)=Check_same(ids, seqs_sub,same,inv)
  for id in seqs:
    if(id in inv):
      if(id in same):
        f = 0
        for id1 in same[id]:
          f = f + sum(map(int, id1.split(read_number_division)[1].split("|")[0].split("_")))
          #f = f + int(id1[len(id1)-1])
        out=out+">"+id.split(read_number_division)[0]+read_number_division+str(f)+"\n"+seqs[id]+"\n"
        ind = ind+1
    else:
      s = seqs[id]
      #if(id.count(read_number_division)>1):
        #id = id.split(read_number_division)
        #id = id[0]+read_number_division+id[len(id)-1]
      out=out+">"+id+"\n"+s+"\n"
      ind = ind+1
    if(ind>100):
      Write_out(out, reduced_sequences)
      (out, ind)=('',0)
  Write_out(out, reduced_sequences)
  return()

def Generate_networks_pre(Sequence_file,tmp_pre,reduced_sequences,edge_lengths,read_number_division,cd_hit_directory):
  Cluster_i(Sequence_file, tmp_pre,edge_lengths,cd_hit_directory)
  (seqs)=Get_seqs_single (Sequence_file)
  Get_sequence_same_single(seqs, tmp_pre,reduced_sequences, tmp_pre,read_number_division)
  return()

def Get_diff(s1,s2,mismatch):
  p=1
  (mm)=Do_counting(s1, s2,mismatch)
  if (mm > mismatch):
    p=0
  return (p,mm)

def Do_counting(s1, s2,mismatch):
  i1=0
  i2=0
  mm=0
  l=min([len(s1),len(s2)])
  for i in range(0,l-2):
    if (s1[i1]==s2[i2]):
      i1=i1+1
      i2=i2+1
    else:
      mm=mm+1
      i2=i2+1
      i1=i1+1
    if (mm>mismatch+1):
      break
  return (mm)

def Generate_networks(Sequence_file, tmp_file1,edge_lengths, tmp_file, file_out,cd_hit_directory):
  Cluster_i(Sequence_file, tmp_file1,edge_lengths,cd_hit_directory)
  (s_sizes, cluster)=Get_cluster_sizes_single (tmp_file1)
  (seqs)=Get_seqs_single (Sequence_file)
  Get_similar_clusters(s_sizes, cluster, seqs,tmp_file+"_coclustered")
  (inv,coclust)=Get_coclustered (tmp_file+"_coclustered")
  Get_cluster_similarities_single(seqs,coclust, cluster,file_out,inv)
  return()

def Get_cluster_similarities_single(seqs,coclust, cluster,file_out,inv):
  fh=open(file_out,"w")
  fh.close()
  ind=0
  total=len(cluster)
  t=0 
  for c in cluster:
    if (c not in inv):
      ind=ind+1
      clust_seqs=[]
      t=0 
      for id in cluster[c]:
        clust_seqs.append((id, seqs[id.split("__")[0]], len(seqs[id.split("__")[0]])))
        t=t+1
      for c1 in coclust[c]:
        for id in cluster[c1]:
          clust_seqs.append((id, seqs[id.split("__")[0]], len(seqs[id.split("__")[0]])))
          t=t+1
      clust_seqs=sorted(clust_seqs, key=itemgetter(2), reverse=True)
      if(len(clust_seqs)>1):
        if(len(clust_seqs)>500):print len(clust_seqs), total, ind
        Get_similarity_single(clust_seqs,file_out,c)
  return()

def Get_similarity_single(clust_seqs,file_out,c):
  done=Tree()
  total = len(clust_seqs)
  out='' 
  ind1 = 0 
  mismatch = 1
  for i1 in range(0,total-1):
    read1 = clust_seqs[i1]
    id1 = read1[0]
    seq1= read1[1]
    l1 = read1[2]
    done[id1]=1
    for i2 in range(i1+1,total):
      if (i1<i2):
        read2 = clust_seqs[i2]
        id2 = read2[0]
        seq2= read2[1]
        l2 = read2[2]
        if (id2 not in done):
          (s1, s2, p1)=Trim_sequences(seq1,seq2,l1,l2)
          if (p1==1):
            if(s1==s2):
              out=out+"0\t"+id1+"\t"+id2+"\t"+c+"\t"+str(l1)+"\t"+str(l2)+"\n"
              ind1 =ind1+1
            else:
              (p,mm)=Get_diff(s1,s2,mismatch)
              if (p==1 and mm<=mismatch):
                out=out+str(mm)+"\t"+id1+"\t"+id2+"\t"+c+"\t"+str(l1)+"\t"+str(l2)+"\n"
                ind1 = ind1+1
        if(ind1>=100):
          Write_out(out, file_out)
          ind1=0
          out=''
  Write_out(out, file_out)
  del out
  return()

def Get_cluster_sizes_single (file):
  (cluster)= Get_clusters(file+".bak.clstr")
  print len(cluster)
  sizes=[] 
  for c in cluster:
    tot=len(cluster[c])
    t=(c, tot)
    if (tot>50):
      sizes.append(t)
  s=sorted(sizes, key=itemgetter(1), reverse=True)
  return(s, cluster)

def Get_vaguely_similar_seqs(s1,s2,mis):
  l1=len(s1)
  l2=len(s2)
  trim=15
  seg_length = (l1-(2*trim))/mis
  p=0
  for i in range(0,mis):
    pos=trim+(i*seg_length)
    seg = s1[pos:pos+seg_length]
    index=s2.find(seg)
    if (index!=-1):
      if(index>pos):
        s2=s2[index-pos:l2]
      else:
        s1=s1[pos-index:l1]
      min_len=min([len(s1),len(s2)])
      s1=s1[0:min_len]
      s2=s2[0:min_len]
      p=1
      break
  return (s1, s2, p)

def Count_diffs (s1, s2, mis):
  i1=0
  i2=0
  mm=0
  p=1
  for i in range(0,len(s2)-1):
    if (s1[i1]==s2[i2]):
      i1=i1+1
      i2=i2+1
    else:
      if (s1[i1+1]==s2[i2]):
        i1=i1+1
        mm=mm+1
      else:
        if (s1[i1]==s2[i2+1]):
          i2=i2+1
          mm=mm+1
        else:
          mm=mm+1
    if (mm>mis):
      p=0
      break
  return (mm,p)

def Get_similar_clusters(s_sizes, cluster, seqs,tmp_file):
  coclust=Tree()
  inv={}
  mis = 5
  out = ''
  indw = 0
  fh1=open (tmp_file, "w")
  fh1.close()
  comp = 5
  for i in range(0,len(s_sizes)):
    print len(s_sizes), i
    clust=s_sizes[i]
    c1=clust[0]
    if (c1 not in inv):
      s1=''
      coclust[c1][c1].value=1
      inv[c1]=c1
      seqs1= []
      ind = 0
      for id in cluster[c1]:
        if (id in seqs):
          s1=seqs[id]
          seqs1.append(s1)
          ind = ind+1
          if(ind>comp):
            break
      for i2 in range(i+1,len(s_sizes)):
        clust2=s_sizes[i2]
        c2=clust2[0]
        ind = 0
        found = 0
        for id in cluster[c2]:
          if (id in seqs):
            ### Get difference with c1 sequences
            ### If sequence is different > X% then move onto next cluster
            ### If sequence is less different than Y% then co-cluster the clusters
            p=0
            ind = ind+1
            for s1 in seqs1:
              (s1, s2,p)=Get_vaguely_similar_seqs(s1, seqs[id], mis)
              if (p==1):
                (mm, q)=Count_diffs(s1, s2,mis)
                if (q==1):
                  out=out+c1+"\t"+c2+"\n"
                  indw = indw+1
                  if (indw >100):
                    Write_out(out, tmp_file)
                    out = ''
                    indw=0
                  found = 1
                  break
            if (ind > comp or found ==1):
              break
  Write_out(out, tmp_file)
  print out
  return()

def Get_coclustered (file):
  inv={}
  coclust = Tree()
  fh=open (file,"r")
  for l in fh:
    l=l.strip()
    l=l.split()
    inv[l[1]]=l[0]
    coclust[l[0]][l[1]].value = 1
  fh.close()
  return(inv, coclust)

def Decon_edges(att_file,file_seqs,file_edges,file_vertex,seq_file):
  fh = open (seq_file , "r")
  freq_id={}
  for header,sequence in fasta_iterator(fh):
    freq_id[header.split("__")[0]]=header
  fh.close()
  inverse,raw=Get_inverse_ids(file_seqs,file_vertex,freq_id)
  edges,edges23 = Tree(), Tree()
  fh1 = open (att_file , "r")
  for l in fh1:
    l=l.strip()
    l1=l.split()
    if (int(l1[0])==1 or int(l1[0])==2):
      id1, id2 = l1[1].split("__")[0],l1[2].split("__")[0]
      id1, id2 = freq_id[id1],freq_id[id2]
      edges[id2][id1].value=1
  fh1.close()
  Print_single_edges(file_edges, inverse, edges,tmp_file1,raw)
  del inverse, edges, edges23
  return()

def Get_inverse_ids(file_seqs,file_vertex,freq_id):
  fh = open (file_seqs, "r")
  inverse,raw = {},{}
  ind =  0
  for l in fh:
    if (l[0]==">"):
      l=l.strip()
      l=l.replace(">","")
      l=l.split("||")
      inverse[freq_id[l[0].split("__")[0]]]=freq_id[l[1].split("__")[0]]
      raw[l[1].split("__")[0]] = 1
      ind = ind+1
  fh.close()
  fh=open(file_vertex,"r")
  for l in fh:
    l=l.strip().split()
    if(l[0].split("__")[0] in raw):
      raw[l[0].split("__")[0]] = l[0]
  fh.close()
  return(inverse,raw)

def Reduce_barcode_forward_reverse_multiplicity(primer_tag_file,Output_trim):
  fh = open(primer_tag_file,"r")
  bcs,seqs,uniq_bc = Tree(),{},{}
  bc_inv = Tree()
  for l in fh:
    if(l[0]!="#"):
      l=l.strip().split()
      if(l[6]=="YES"):
        bc1,bc2,consensus,id,n_reads = l[3],l[4],l[7],l[0], int(l[2])
        seqs[id]  = [consensus, n_reads]
        bcs[bc1][bc2][id].value = 1
        bcs[bc2][bc1][id].value = 1
        bc_inv[consensus][bc1][bc2][id].value = 1
        bc_inv[consensus][bc2][bc1][id].value = 1
        uniq_bc[bc1]=1
        uniq_bc[bc2]=1
  fh.close()
  same = Tree()
  for s in bc_inv:
    for bc1 in bc_inv[s]:
      ids = []
      for bc2 in bc_inv[s][bc1]:
        for id in bc_inv[s][bc1][bc2]:
          ids.append(id)
        if(len(ids)>0):
          for i in range(0,len(ids)):
            for j in range(i,len(ids)):
              if(i<j):
                same[ids[i]][ids[j]].value = 1 
                same[ids[j]][ids[i]].value = 1 
  sub_same, sub_inv=Deconvolute_same_array (same)
  del same, bcs
  primer_tag_file_post = primer_tag_file.replace(".txt","_REDUCED.txt")
  for f in [primer_tag_file_post,Output_trim]:
    fh=open(f,"w")
    fh.close()
  out, ind= '#Id\tsequence\tmultiplicity\n', 0
  all_seqs =Tree()
  for i in seqs:
    seq=seqs[i][0]
    seq_sub = {}
    if i in sub_same:
      for j in sub_same[i]:
        f = seqs[j][1]
        s1 = seqs[j][0]
        if(s1 in seq_sub):seq_sub[s1] = seq_sub[s1]+f
        else: seq_sub[s1]=f
    else:
      seq_sub[seqs[i][0]] = seqs[i][1]
    for s in seq_sub:
      out=out+i+"\t"+seq+"\t"+str(seq_sub[s])+"\n"
      all_seqs[seq][i].value = 1
      ind = ind+1
      if(ind>500):
        Write_out(out, primer_tag_file_post)
        out, ind = '',0
  Write_out(out, primer_tag_file_post)
  out, ind = '',0
  for s in all_seqs:
    for id in all_seqs[s]:
      break
    out=out+">"+id.split("__")[0]+"__"+str(len(all_seqs[s]))+"\n"+s+"\n"
    ind = ind+1
    if(ind>500):
      Write_out(out,Output_trim)
      out, ind = '',0
  Write_out(out,Output_trim)
  return()
    
def Decon_identical (seq_file, att_file, file_vertex, file_seqs, tmp_file,read_number_division):
  fh = open (seq_file , "r")
  all,seqs,freq_id ={},{},{}
  for header,sequence in fasta_iterator(fh):
    seqs[header.split("__")[0]]=sequence
    all[header.split("__")[0]]=sequence
    freq_id[header.split("__")[0]]=header
  fh.close()
  fh1 = open (att_file , "r")
  same1 = Tree()
  for l in fh1:
    l=l.strip()
    l1=l.split()
    cluster = l1[3]
    if (l1[0]=="0"):
      same1[cluster][l1[2]][l1[1]].value = 1
      same1[cluster][l1[1]][l1[2]].value =1 
  fh1.close()
  same,inverse, inv, out, ind, length = {},{},{},'',0,{}
  fh=open (file_seqs, "w")
  fh1.close()
  j=header
  for c in same1:
    sub_same = Tree()
    for id1 in same1[c]:
      for id2 in same1[c][id1]:
        sub_same[id1.split("__")[0]][id2.split("__")[0]].value=1
        sub_same[id2.split("__")[0]][id1.split("__")[0]].value=1
    (sub_same, sub_inv)=Deconvolute_same_array (sub_same)
    for i in sub_same:
      s=seqs[i.split("__")[0]]
      total = 0
      mins=s
      for j in sub_same[i]:
        j1 = freq_id[j.split("__")[0]]
        freq = map(int, j1.split(read_number_division)[1].split("|")[0].split("_"))
        if(total==0):total = freq
        else:
          total = [a + b for a, b in zip(freq,total)]
          #total = map(add,freq,total)
        inverse[j.split("__")[0]]=i
        out = out+">"+j+"||"+i+"\n"+seqs[j.split("__")[0]]+"\n"
        s=seqs[j.split("__")[0]]
        if(len(s)<len(mins)):
          mins=s
        ind = ind+1
        if(ind>100):
          Write_out(out, file_seqs)
          out = ""
          ind = 0
      same[i]=total
      length[i]=mins
  info = ''
  if(len(freq_id[j.split("__")[0]].split(read_number_division)[1].split("|"))>=2):
    info = "|"+freq_id[j.split("__")[0]].split(read_number_division)[1].split("|")[1]
  Write_out(out, file_seqs)
  del seqs
  Print_vertices(all, inverse, same, file_vertex, length,read_number_division,info,freq_id)
  return()

def Print_single_edges(file_edges, inverse, edges,tmp_file,raw):
  fh = open (tmp_file, "w")
  fh.close()
  edge,ind = '',0
  for id1 in edges:
    ida=id1 
    if (id1 in inverse):
      ida = inverse[id1]
      if(ida.split("__")[0] in raw):
        ida = raw[ida.split("__")[0]]
    for id2 in edges[id1]:
      idb = id2
      if (id2 in inverse):
        idb = inverse[id2]
        if(idb.split("__")[0] in raw):
          idb = raw[idb.split("__")[0]]
      if (ida!=idb):
        edge=edge+ida+"\t"+idb+"\t"+str(1)+"\t"+id1+"\t"+id2+"\n"
        ind = ind+1
        if (ind>300):
          Write_out(edge, tmp_file)
          edge = ''
          ind = 0
  Write_out(edge, tmp_file)
  del edges
  return()

def Print_vertices(all, inverse, same, file_vertex,length,read_number_division,info,freq_id):
  out =''
  total = 0
  fh = open (file_vertex, "w")
  fh.close()
  ind =0
  for id in all:
    ind = ind+1
    if (id in inverse):
      if (id in same):
        id1 = id.split(read_number_division)[0]+read_number_division+"_".join(map(str,same[id]))+info
        if(id not in length):
          out = out +id1+"\t"+str(same[id])+"\t"+all[id]+"\n"
        else:
          out = out +id1+"\t"+str(sum(same[id]))+"\t"+length[id]+"\n"
    else:
      if(id in freq_id):
        id1 =freq_id[id]
      else:print id
      freq = sum(map(int, id1.split(read_number_division)[1].split("|")[0].split("_")))
      out = out +id1+"\t"+str(freq)+"\t"+all[id]+"\n"
    if (ind>300):
      Write_out(out, file_vertex)
      ind =0
      out = ''
  Write_out(out, file_vertex)
  return()

def Deconvolute_same_array (tree):
  decon = Tree()
  inv = {}
  index = 0
  for i in tree:
    if (i not in inv):
      index = index+1
      decon[i][i].value = 1
      inv[i]=i
      array = []
      for j in tree[i]:
        decon[i][j].value = 1
        inv[j]=i
        array.append(j)
      array_new=[]
      found = 0
      if (len(array)>0):
        found = 1
      while (found == 1):
        array_new =[]
        for k in array:
          for l in tree[k]:
            if (l not in inv):
              array_new.append(l)
              inv[l]=i
              decon[i][l].value = 1
        array = array_new
        if (len(array_new)==0):
          found = 0
  return (decon, inv)

def Reduce_edges(file_in, file_out):
  done,out, ind = Tree(),'',0
  fh = open (file_in, "r")
  fh1 = open (file_out, "w")
  fh1.close()
  for l in fh:
    l=l.strip()
    l=l.split()
    if (l[0] not in done[l[1]] and l[1] not in done[l[0]]):
      if (l[0]!=l[1]):
        out = out+l[0]+"\t"+l[1]+"\t"+l[2]+"\n"
        done[l[0]][l[1]].value =1
        ind = ind+1
        if (ind>300):
          ind =0
          Write_out(out, file_out)
          out=''
  Write_out(out, file_out)
  del done
  return()

def Deconvolute_edges(seq_file,  att_file, file_vertex, file_seqs, tmp_file0, file_edges,read_number_division):
  Decon_identical (seq_file, att_file, file_vertex, file_seqs, tmp_file0,read_number_division)
  Decon_edges(att_file,file_seqs,file_edges,file_vertex,seq_file)
  Reduce_edges (tmp_file1, file_edges)
  return()

def Get_network_clusters(file_vertex, file_edges, cluster_file):
  (G,scale)=Read_graphical_inputs(file_vertex, file_edges)
  Output_cluster_file (G,cluster_file)
  return()

def Read_graphical_inputs(file_vertex, file_edges):
  fh = open (file_vertex, "r")
  freq = {}
  size=[]
  ind = 0
  G=nx.Graph()
  G.rtt={}
  for l in fh:
    l=l.strip()
    l=l.split()
    f = int(l[1])#sum(map(int, l[1].split(read_number_division)[1].split("|")[0].split("_")))
    freq[l[0]]=f
    size.append(f)
    G.add_node(l[0])
    G.rtt[l[0]] = int(f)
    ind = ind+1
  scale1 = max(size)
  fh.close()
  fh = open (file_edges, "r")
  for l in fh:
    l=l.strip()
    l=l.split()
    if (l[0] in freq and l[1] in freq):
      G.add_edge(l[0],l[1])
  fh.close()
  return(G,scale1)

def Output_cluster_file (G, cluster_file):
  con= nx.connected_components(G)
  ind = 0
  ind1 = 0
  ind2 = 0
  fh = open(cluster_file, "w")
  fh.close()
  out = '# Connected_components\n'
  max_f,t,nvertmax = 0,0,0
  for i in con:
    ind = ind+1
    tc = 0
    nvert = 0
    for j in i:
      ind1=ind1+1
      ind2 = ind2+1
      out = out+str(ind1)+"\t"+str(ind)+"\t"+j+"\t"+str(G.rtt[j])+"\n"
      tc,t = tc+G.rtt[j],t+G.rtt[j]
      nvert = nvert+1
      if (ind2>100):
        Write_out(out, cluster_file)
        out = ""
        ind2=0
    if(tc>max_f):max_f,nvertmax = tc,nvert
  Write_out(out, cluster_file)
  print file_vertex,"Maximum cluster:",max_f*100.0/t,"%", nvertmax , "vertices"
  return()

def Get_network_input(file_cluster, outfile,edge_file, checked_edges):
  fh=open(outfile,"w")
  fh.close()
  fh=open(file_cluster, "r")
  cluster=Tree()
  ids = {}
  for l in fh:
    if(l[0]!="#"):
      l=l.strip().split()
      cluster[l[1]][l[1]+"\t"+l[2]+"\t"+l[3]].value=1
      ids[l[2]]=1
  fh.close()
  (out, ind)=('',0)
  for c in cluster:
    if(len(cluster[c])>1):
      for l in cluster[c]:
        out=out+l+"\n"
        ind = ind+1
        if(int(l.split("\t")[2])>1000):print l
    else:
      for l in cluster[c]:
        if(int(l.split("\t")[2])>1):
          out=out+l+"\n"
          ind = ind+1
    if(ind>300):
      Write_out(out, outfile)
      (out, ind)=('',0)
  Write_out(out, outfile)
  fh=open(checked_edges,"w")
  fh.close()
  fh=open(edge_file,"r")
  (out, ind)=('',0)
  for l in fh:
    l=l.strip().split()
    if(l[0] in ids and l[1] in ids):
      out=out+l[0]+"\t"+l[1]+"\t"+l[2]+"\n"
      ind = ind+1
      if(ind>300):
        Write_out(out, checked_edges)
        (out, ind)=('',0)
  Write_out(out, checked_edges)
  fh.close()
  return()

def Reduce_identical_sequences(Reduced_file,file_vertex,read_number_division):
  fh=open(Reduced_file,"w")
  fh.close()
  fh=open(file_vertex,"r")
  (ind, out) = (0,'')
  for l in fh:
    l=l.strip().split()
    if(l[0].count("|") >=1):out=out+">"+l[0]+"\n"+l[2]+"\n"
    else:out=out+">"+l[0].split(read_number_division)[0]+read_number_division+l[1]+"\n"+l[2]+"\n"
    ind = ind+1
    if(ind>500):
      Write_out(out, Reduced_file)
      (ind, out) = (0,'')
  fh.close()
  Write_out(out, Reduced_file)
  return()


class Tree(defaultdict):
  def __init__(self, value=None):
    super(Tree, self).__init__(Tree)
    self.value = value

def Write_out(out, file):
  fh=open(file, "a")
  fh.write(out)
  fh.close()
  return ()

def Get_locations(ref_locations):
  fh=open(ref_locations,"r")
  locations = {}
  for l in fh:
    l=l.strip().split()
    locations[l[0]] = l[1]
  fh.close()
  return(locations)

###########################
if(len(sys.argv)<10):
  print "\nUsage:\npython Read_processing_and_quality_PUBLIC.py <output directory> <sample ID> <barcode group> <chain type (IGH/TCR)> <paired type> <species> <source files> <thresold length> <primer file> <method> <optional: reverse primer group>"
  print "\nNote:\n\tPlease check line XXX location of directory with reference libraries\n"
  #quit()

dir = sys.argv[1]
id = sys.argv[2]
barcode_group = sys.argv[3]
gene = sys.argv[4]
paired = sys.argv[5]
species = sys.argv[6]
source = sys.argv[7]
length=sys.argv[8]
primer_file = sys.argv[9]
method = sys.argv[10]
command_source = sys.argv[11]
command_source = command_source.split(",")
ORF_filter = "TRUE"
if(len(sys.argv)>13):
  reverse_primer_group = sys.argv[13]
if(len(sys.argv)>14):
  ORF_filter = sys.argv[14]
else:reverse_primer_group = "OTHER"
print "Reverse primer group: ",reverse_primer_group
other = sys.argv[12]
###########################
pwd = "/well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/"#commands.getoutput("pwd")
ref_locations = pwd+"/Locations_of_called_programmes.txt"
########################### Files for QC and filtering
Seq_file1 = dir+"FASTQ_FILES/Sequences_"+id+"_1.fasta"
Seq_file2 = dir+"FASTQ_FILES/Sequences_"+id+"_2.fasta"
#Tmp_file = dir+"ORIENTATED_SEQUENCES/TMP/Untrimmed_"+id+".fasta"
if(method=="Multiplex_FORWARD_BARCODE_GROUPED" or method=="Multiplex_REVERSE_BARCODE_GROUPED"): 
  Tmp_file = dir+"FASTQ_FILES/Seqs_barcode_split_"+barcode_group+"_"+other.split(",")[0]+".fasta" 
else: Tmp_file = dir+"FASTQ_FILES/Seqs_barcode_split_"+barcode_group+".fasta"
Trim1=dir+"ORIENTATED_SEQUENCES/TMP/Trimmed_orientated_all_"+id+".fasta"
Trim2=dir+"ORIENTATED_SEQUENCES/TMP/Filtered_J_"+id+".fasta"
Trim3=dir+"ORIENTATED_SEQUENCES/TMP/Filtered_reduced_"+id+".fasta"
Fail_file = dir+"FASTQ_FILES/Fail_filtered_"+id+".fasta"
primer_tag_file = dir+"ORIENTATED_SEQUENCES/TMP/Barcode_filtering_information_"+id+".txt"
primer_tag_file_count = dir+"ORIENTATED_SEQUENCES/TMP/All_barcodes_"+id+".txt"
Filtered_out1=dir+"ORIENTATED_SEQUENCES/Filtered_ORFs_sequences_all_"+id+".fasta"
nn_orf_filtered = dir+"ORIENTATED_SEQUENCES/Nucleotide_ORF_filtered_all_"+id+".fasta"
tmp_file_orf = dir+"ORIENTATED_SEQUENCES/TMP/Blast_matching_"+id
filtering_report = dir+"ORIENTATED_SEQUENCES/Filtering_report_"+id+".txt"
########################## Files for clustering
att_file = dir+"ORIENTATED_SEQUENCES/NETWORKS/Vertex_relations_"+id+".txt"
file_vertex = dir+"ORIENTATED_SEQUENCES/NETWORKS/Att_"+id+".txt"
file_edges = dir+"ORIENTATED_SEQUENCES/NETWORKS/Edges_"+id+".txt"
cluster_file = dir+"ORIENTATED_SEQUENCES/NETWORKS/Cluster_identities_"+id+".txt"
Reduced_file = dir+"ORIENTATED_SEQUENCES/NETWORKS/Fully_reduced_"+id+".fasta"
checked_edges = dir+"ORIENTATED_SEQUENCES/NETWORKS/Checked_edges_"+id+".txt"
plot_ids_file = dir+"ORIENTATED_SEQUENCES/NETWORKS/Plot_ids_"+id+".txt"
file_seqs = dir+"ORIENTATED_SEQUENCES/NETWORKS/Sequences_"+id+".txt"
tmp_file0 = dir+"ORIENTATED_SEQUENCES/NETWORKS/Decon_0_"+id+".txt"
tmp_reduced_sequences = dir+"ORIENTATED_SEQUENCES/NETWORKS/Primary_reduced_sequences"+id+".fasta"
tmp_pre = dir+"ORIENTATED_SEQUENCES/NETWORKS/Pre_tmp_"+id
tmp_file1 = dir+"ORIENTATED_SEQUENCES/NETWORKS/NN_Tmp_cluster_"+id+".1"
edge_lengths=0.85
tmp_file = dir+"ORIENTATED_SEQUENCES/NETWORKS/NN_Tmp_cluster_"+id+"."
read_number_division = "__"
########################## Reference files
locations = Get_locations(ref_locations)
gene_sub = gene
if(gene in ["TRB"]):gene_sub = "TCR"
quasr_location = locations["quasr_location"] 
ref_dir = locations["ref_library"]
refv = ref_dir+"Reference_nn_"+species+"_"+gene+"V.fasta"
refj = ref_dir+"Reference_nn_"+species+"_"+gene+"J.fasta"
refvp= ref_dir+"Reference_protein_"+species+"_"+gene+"V.fasta"
refjp= ref_dir+"Reference_protein_"+species+"_"+gene+"J.fasta"
ref_const = ref_dir+"Reference_nn_"+species+"_"+gene_sub+"_constant_exon1.fasta"
cd_hit_directory = locations["cd_hit_directory"]
blastall = locations["blastall"]
readjoin_location = locations["read_join"]
codon_table = locations["codon_table"]
primer_length = 20
######################### Commands
if(command_source.count("1")!=0):
  Intialise_files(dir)
  QC_samples(dir,gene,id,source,length,species,barcode_group,quasr_location,readjoin_location)
  Split_reads_and_join(Seq_file1, Seq_file2, Tmp_file,gene,paired,id,method,primer_file,barcode_group,sys.argv[12],dir,ref_const)
  
### Tip: it is good to check all the fasta files in the FASTQ_FILES directory have been made correctly at this point (with non-zero number of lines)

print gene
### Filtering and processing reads
if(command_source.count("2")!=0):
  if(gene.count("IG")!=0):
    Trim_sequences_BCR_TCR(Tmp_file,Fail_file,Trim1, gene,paired,species,primer_file,primer_tag_file,tmp_file,primer_tag_file_count,id,ref_const,reverse_primer_group,cd_hit_directory,other,barcode_group)
    Filter_IgJ_genes(Trim1, Trim2, refj,primer_file,ref_const,primer_tag_file_count,blastall,refv,primer_length)
    Reduce_sequences(Trim2, Trim3,primer_file)
    ORF_calculation_single (Trim3, Filtered_out1,nn_orf_filtered,dir,gene,refv, refj,refvp,refjp,tmp_file_orf,codon_table,ORF_filter)
    Get_read_report(Seq_file1, Seq_file2, Tmp_file, Trim1, nn_orf_filtered,filtering_report,id, species, gene, dir,primer_tag_file_count,primer_file,method,barcode_group,Trim3,other)
  if(gene.count("TR")!=0 or gene.count("TCR")!=0): #### NOT RIGHT FOR TCR THAT ARE NOT BARCODED
    Trim_sequences_BCR_TCR(Tmp_file,Fail_file,Trim1, gene,paired,species,primer_file,primer_tag_file,tmp_file,primer_tag_file_count,id,ref_const,reverse_primer_group,cd_hit_directory,other,barcode_group)
    Filter_IgJ_genes(Trim1, Trim2, refj,primer_file,ref_const,primer_tag_file_count,blastall,refv,primer_length)
    Reduce_sequences(Trim2, Trim3,primer_file)
    ORF_calculation_single (Trim3, Filtered_out1,nn_orf_filtered,dir,gene,refv, refj,refvp,refjp,tmp_file_orf,codon_table,ORF_filter)
    Get_read_report(Seq_file1, Seq_file2, Tmp_file, Trim1, nn_orf_filtered,filtering_report,id, species, gene, dir,primer_tag_file_count,primer_file,method,barcode_group,Trim3,other)

### Clustering reads
if(command_source.count("3")!=0):
          #### Clustering: 1. Reduce large datasets into unique sequences 
          #### Clustering: 2. Cluster related sequences together
  Generate_networks(nn_orf_filtered, tmp_file1,edge_lengths, tmp_file, att_file,cd_hit_directory)
  Deconvolute_edges(nn_orf_filtered,  att_file, file_vertex, file_seqs, tmp_file0, file_edges,read_number_division)
  Get_network_clusters(file_vertex, file_edges, cluster_file)
  Reduce_identical_sequences(Reduced_file,file_vertex,read_number_division)
  Get_network_input(cluster_file,plot_ids_file,file_edges, checked_edges)


