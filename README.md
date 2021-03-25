# BCR_TCR_PROCESSING_PIPELINE
BCR/TCR processing pipeline for NGS data

Please see the manual for full details. A summary of the pipeline is described below:
Page 2 of 11
 
**Stage 1:**
1. QC sequences

**Stage 2:**
1. Join forward and reverse reads (merging)
2. Split sequences according to sample barcode
3. Identify RNA barcode and collapse/error correct based on groups of sequences
sharing same barcode
4. Check isotype against reference
5. Check matches to IGHV/J reference sequences
6. Check open reading frame present

**Stage 3:**
1. Network generation: Here, each vertex represents a different sequence, and the number of identical BCR sequences defines the vertex size. Edges are created between vertices that differ by one nucleotide. Clusters are groups of interconnected vertices (1, 2). The program described here calculates edges between unique sequences and determines vertex sizes, creating output files in formats that can directly be used in network analysis programs such as networkx (python) or igraph (R or python).

**Stage 4:**
1. Sequence annotation
2. Network Analysis
3. Generation of broad repertoire statistics

**Stage 5 (optional):**
1. Run the optional R script to analyse and visualise summary metrics from the results of Stages 1-4. This will also concatenate files for all samples.
2. Check the percentage of reads which pass Open-Read-Frame filtering.
a. If this is abnormally low (potentially due to large clonal expansion) consider running the analysis with the ORF column in the sample sheet set to
anything other than TRUE. To prevent reads being removed.

**Stage 6:**
1. Concatenate filtered fastq files into a smaller number of multi-individual large files.
2. Upload large files to IMGT for annotation.

**Stage 7:**
1. Results of IMGT analysis are used in Isotyper specific Analysis
