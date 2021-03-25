#!/usr/bin/python
#import math
import sys
import os
import commands

def Get_info(file):
  (id,  sample,  info,  gene,  directory,pair,pair_final,dir, platform,spec,primer,other,reverse_primer_group)=([],[],[],[],[],[],[],[],[],[],[],[],[])
  fh=open(file,"r")
  ind = 0
  for l in fh:
    ind = ind+1
    if(l[0]!="#" and len(l)>3):
     # l=l.replace('/lustre/scratch108/viruses/rbr1/','/lustre/scratch118/infgen/team146/rbr1/')
      l=l.replace("/lustre/scratch118/infgen/team146/rbr1/MISEQ/","/well/immune-rep/users/mfj169/MISEQ/")
      l=l.strip().split()
      id.append(l[0])
      sample.append(l[1])
      info.append(l[2])
      gene.append(l[3])
      directory.append(l[4])
      pair.append(l[5])
      pair_final.append(l[6])
      dir.append(l[7])
      platform.append(l[8])
      spec.append(l[9])
      if(len(l)>=11):primer.append(l[10])
      else:primer.append("LIBRARY/FR1_primers.txt")
      if(len(l)>=14):reverse_primer_group.append(l[13])
      else:reverse_primer_group.append("STANDARD")
      if(len(l)>=12):
        lis = []
        for i in range(11,len(l)):
          lis.append(l[i])
        other.append(",".join(lis))
      else:other.append('')
  fh.close()
  return(id,  sample,  info,  gene,  directory, pair,pair_final,dir, platform,spec,primer,other,reverse_primer_group)

args=sys.argv
queue = "short.qc"
#queue = "long.qc"
python = "/usr/bin/python"
if (len(args)<5):
  queue = '-q normal'
  command = 'python Processing_sequences_large_scale.py [sample file list] [commands (comma separated list)] [bsub command: Y/N] [print commands: Y/N] [run commands: Y/N]'
  print "SEQUENCE ANALYSIS PIPELINE: Creates networks from MiSeq data"
  print "USAGE:"
  print command,"\n"
  os.system("cat Command_outline.txt")
  print "\n"
else:
  file = args[1]
  command = args[2]
  bsub_command = args[3]
  print_command = args[4]
  run_command = args[5]
  command=command.split(",")
  (ids,  samples,  infos,  gene,  source,pairs,pairs_final,dirs, platform,spec,primer,others,reverse_primer_group)=Get_info(file)
  wkg_dir = commands.getoutput("pwd")+"/"
  if(bsub_command not in ["Y","N"]):print 'python Processing_sequences.py [sample file list] [commands (comma separated list)] [bsub command: Y/N] [print commands: Y/N] [run commands: Y/N]\n\tError: bsub command must be: Y or N'
  if(bsub_command not in ["Y","N"]):print 'python Processing_sequences.py [sample file list] [commands (comma separated list)] [bsub command: Y/N] [print commands: Y/N] [run commands: Y/N] \n\tError: print command must be: Y or N'
  if(bsub_command not in ["Y","N"]):print 'python Processing_sequences.py [sample file list] [commands (comma separated list)] [bsub command: Y/N] [print commands: Y/N] [run commands: Y/N] \n\tError: run command must be: Y or N'
  idss,dirss='',''
  commands = []
  bsubs = []
  for i in range(0,len(samples)):
    info,sample, gene_types, pair,pair_final, dir,platforms, primers,other=infos[i], samples[i], gene[i], pairs[i],pairs_final[i], dirs[i], platform[i], primer[i],others[i]
    bsub = ''
    id,sources,species=ids[i],source[i],spec[i]
    if(bsub_command=="Y"):
      bsub = " | xargs -i echo qsub  -P immune-rep.prjc -q "+queue+" -b y -o "+wkg_dir+"out_STANDARD_"+id+"  -e "+wkg_dir+"error_log_"+id+" -N job_name \"{}\"  | sh"
    bsubs.append(bsub)
    if( "1" in command):
      command1 = python+" "+wkg_dir+"Read_processing_and_quality.py "+dir+" "+id+" "+sample+" "+gene_types+" "+pair+" "+species+" "+sources +" "+str(200)+" "+primers+" "+platforms+" 1 "+other+" "+reverse_primer_group[i]
      commands.append(command1)
    if( "2" in command):
      command1 = python+" "+wkg_dir+"Read_processing_and_quality.py "+dir+" "+id+" "+sample+" "+gene_types+" "+pair+" "+species+" "+sources +" "+str(200)+" "+primers+" "+platforms+" 2 "+other+" "+reverse_primer_group[i]
      commands.append(command1)
    if( "3" in command):
      command1 = python+" "+wkg_dir+"Read_processing_and_quality.py "+dir+" "+id+" "+sample+" "+gene_types+" "+pair+" "+species+" "+sources +" "+str(200)+" "+primers+" "+platforms+" 3 "+other+" "+reverse_primer_group[i]
      commands.append(command1)
    if( "4" in command):
      command2 = python+" "+wkg_dir+"Generate_repertoire_statistics.py "+dir+"ORIENTATED_SEQUENCES/ANNOTATIONS/ "+id+" "+dir+"ORIENTATED_SEQUENCES/NETWORKS/Fully_reduced_"+id+".fasta "+dir+"ORIENTATED_SEQUENCES/Filtered_ORFs_sequences_all_"+id+".fasta "+gene_types+" "+species+" "+dir+"ORIENTATED_SEQUENCES/NETWORKS/Cluster_identities_"+id+".txt ANNOTATE,STATISTICS "+reverse_primer_group[i]
      commands.append(command2)
    if("ISO1" in command):
      command1 = "python "+wkg_dir+"IsoTyper_1.0.py "+id+" "+id+" "+dir+" "+species+" "+reverse_primer_group[i]
      commands.append(command1)
    if("ISO2" in command):
      command1 = "python "+wkg_dir+"Per_isotype_cluster_analyses.py "+id+" "+id+" "+dir+" "+species+" "+reverse_primer_group[i]
      commands.append(command1)
    if("NONISO1" in command):
      command1 = "python "+wkg_dir+"Non_isotyper_1.0.py "+id+" "+id+" "+dir+" "+species+" "
      commands.append(command1)
    if("TCRISO1" in command):
      command1 = "python "+wkg_dir+"TCRoTyper_1.0.py "+id+" "+id+" "+dir+" "+species+" "+reverse_primer_group[i]
      commands.append(command1)
    if("CSR" in command):
      command14 = "python "+wkg_dir+"Class_switch_recombination_analysis.py "+dir+"ORIENTATED_SEQUENCES/ANNOTATIONS/ "+id+" "+dir+"ORIENTATED_SEQUENCES/NETWORKS/Fully_reduced_"+id+".fasta "+dir+"ORIENTATED_SEQUENCES/Filtered_ORFs_sequences_all_"+id+".fasta "+gene_types+" "+species+" "+dir+"ORIENTATED_SEQUENCES/NETWORKS/Cluster_identities_"+id+".txt 1"
      commands.append(command14)
    if("SUBSAMPLE" in command):
      command16 = "python "+wkg_dir+"Subsampling_networks.py "+dir+"ORIENTATED_SEQUENCES/ANNOTATIONS/ "+id+" "+dir+"ORIENTATED_SEQUENCES/NETWORKS/Fully_reduced_"+id+".fasta "+dir+"ORIENTATED_SEQUENCES/NETWORKS/Edges_"+id+".txt"
      commands.append(command16)
  for i in range(0,len(commands)):
    comm, bsub = commands[i], bsubs[i]
    if(bsub!=''):comm = "echo \'"+comm+"\' "+bsub
    if(print_command=="Y"):print comm, "\n"
    if(run_command=="Y"): os.system(comm)


