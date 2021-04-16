## Script transcritomics
## A. Balard
## April 2021
## stored in Harriet (Heitlinger team server)
TOR_PATH=/SAN/Alices_sandpit/Torelli_transcriptomics

# Kidney cells from Myodes glareolus (BVK), lung cells from Apodemus agrarius (AAL), and kidney cell line from Microtus arvalis (FMN)were either treated with IFNg or not. 
# RNAseq was performed for the 6 samples using Illumina-HiSeq2500/4000. 42.5 million reads were obtained on average per sample (min: 38 million; max: 56.5 million). 

#### SOURCE: https://informatics.fas.harvard.edu/best-practices-for-de-novo-transcriptome-assembly-with-trinity.html

# Part I. CLEAN READS

## Step 1. fastqc, quality metrics reads

# raw reads are in $TOR_PATH/big_data/transcriptomicsTorellidata/1.rawreads
# fastqc reports are in html format in $TOR_PATH/big_data/transcriptomicsTorellidata/2.fastqc_rawreads
# Results: some K-mers beginning/end of reads; quite a bunch of duplicated reads; weird base content at extremities, etc. 

## Step 2. remove erroneous k-mers BEFORE assembly using rcorrector (Song, L. et al. 2015)

# tools: rcorrector (/user/local/bin/ in harriet) on symlinks
cd $TOR_PATH/big_data/transcriptomicsTorellidata/3.filtrated_RNA_reads
perl /usr/local/bin/rcorrector/run_rcorrector.pl -t 4 -1 $TOR_PATH/1.rawreads/AAL-ctrl_1.fq -2 rawRNAseqReads/AAL-ctrl_2.fq
## done for all 6 samples; output format: AAL-ctrl_1.cor.fq
#report in: $TOR_PATH/big_data/transcriptomicsTorellidata/3.filtrated_RNA_reads/rcorrector_report.txt

## Step 3. discard read pairs for which one of the reads is deemed unfixable

#reads are often riddled with Ns, or represent other low complexity sequences. Use https://github.com/harvardinformatics/TranscriptomeAssemblyTools/blob/master/FilterUncorrectabledPEfastq.py
#tool: alice@harriet:/usr/local/bin/FilterUncorrectabledPEfastq.py
cd $TOR_PATH/big_data/transcriptomicsTorellidata/3.filtrated_RNA_reads
sh /usr/local/bin/FilterUncorrectabledPEfastq.py AAL-ctrl_1.cor.fq AAL-ctrl_2.cor.fq 
# done for all 6 samples
# output format: fixed_AAL-ctrl_1.cor.fq

## Step 4. trim adapter and low quality bases from fastq files
#tool: /usr/local/bin/trim_galore_zip/trim_galore (http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)
cd $TOR_PATH/big_data/transcriptomicsTorellidata/4.trimmed_reads
1=$TOR_PATH/big_data/transcriptomicsTorellidata/3.filtrated_RNA_reads/fixed_AAL-ctrl_1.cor.fq 
2=$TOR_PATH/big_data/transcriptomicsTorellidata/3.filtrated_RNA_reads/fixed_AAL-ctrl_2.cor.fq
/usr/local/bin/trim_galore_zip/trim_galore --paired --retain_unpaired --phred33 --output_dir trimmed_reads --length 36 -q 5 --stringency 1 -e 0.1 $1 $2
# done for all 6 samples

## Step 5. Run fastqc on your processed reads that pass qc and filtering from the above steps
# all in $TOR_PATH/big_data/transcriptomicsTorellidata/5.fastqcAfterTrimming

############################################
# PART II. Trasncriptome assembly
## Step 1. Run Trinity DE NOVO
## (A) one transcriptome per SAMPLE 
## (B) one transcriptome by SPECIES (merge Ctrl and IFNg)
## (C) genome-guided

TRINITY_HOME=/usr/local/bin/trinityrnaseq-Trinity-v2.6.6/

# (A) one transcriptome per SAMPLE
# $1 = comma-separated list of R1 files
# $2 = comma-separated list of R2 files
# $3 = name of output directory Trinity will create to store results. This must include Trinity in the name, otherwise the job will terminate
$TRINITY_HOME/Trinity --seqType fq --SS_lib_type RF --max_memory 20G --min_kmer_cov 1 --CPU 4 --left $TOR_PATH/big_data/transcriptomicsTorellidata/4.trimmed_reads/fixed_AAL-ctrl_1.cor_val_1.fq --right $TOR_PATH/big_data/transcriptomicsTorellidata/4.trimmed_reads/fixed_AAL-ctrl_2.cor_val_2.fq --output $TOR_PATH/big_data/transcriptomicsTorellidata/6.TrinityAssemblies/A.oneTranscriptomePerSample/AAL-ctrl
# done 6 times

## (B) one transcriptome by SPECIES (merge Ctrl and IFNg)
$TRINITY_HOME/Trinity --seqType fq --SS_lib_type RF --max_memory 20G --min_kmer_cov 1 --CPU 4 --left $TOR_PATH/big_data/transcriptomicsTorellidata/4.trimmed_reads/merged_by_species/AAL-trimmed-merged-1.fq --right $TOR_PATH/big_data/transcriptomicsTorellidata/4.trimmed_reads/merged_by_species/AAL-trimmed-merged-2.fq --output $TOR_PATH/big_data/transcriptomicsTorellidata/6.TrinityAssemblies/B.oneTranscriptomePerSpecies/AAL-merged_trinity
## RESULT IS $TOR_PATH/big_data/transcriptomicsTorellidata/6.TrinityAssemblies/B.oneTranscriptomePerSpecies/AAL-merged_trinity/AAL-Trinity-transcriptome.fasta
# done 3 times

## (C) genome guided for myodes and microtus
# ?? FMN (microtus arvalis) and BVK (myodes glareolus) guided by microtus ochrogaster genome (https://www.ncbi.nlm.nih.gov/genome/10848)


## Step 2. Quality check
## Transcriptomes quality check

### 1. Assembly metrics were calculating using TrinityStats.pl script 
#### for oneTranscriptomePerSpecies
cd $TOR_PATH/big_data/transcriptomicsTorellidata/7.transcriptomeQualityCheck/B.forOneTranscriptomePerSpecies
$TRINITY_HOME/util/TrinityStats.pl $TOR_PATH/big_data/transcriptomicsTorellidata/6.TrinityAssemblies/B.oneTranscriptomePerSpecies/AAL-merged_trinity/AAL-Trinity-transcriptome.fasta > AAL_stats.txt
# done 3 times

#### for genome guided 

#### TO DO

### 2. Quantification of reads supports was done using bowtie2, samtools and trinity (ref) were used to align the trimmed reads back to their corresponding transcriptome
# 2.1, build a bowtie2 index for your assembly.
bowtie2-build $TOR_PATH/big_data/transcriptomicsTorellidata/6.TrinityAssemblies/B.oneTranscriptomePerSpecies/AAL-merged_trinity/AAL-Trinity-transcriptome.fasta $TOR_PATH/big_data/transcriptomicsTorellidata/6.TrinityAssemblies/B.oneTranscriptomePerSpecies/AAL-merged_trinity/AAL-Trinity-transcriptome

# 2.2. map the reads.
bowtie2 -p 16 --local --no-unal -x $TOR_PATH/big_data/transcriptomicsTorellidata/6.TrinityAssemblies/B.oneTranscriptomePerSpecies/AAL-merged_trinity/AAL-Trinity-transcriptome -q -1 $TOR_PATH/big_data/transcriptomicsTorellidata/4.trimmed_reads/merged_by_species/AAL-trimmed-merged-1.fq -2 $TOR_PATH/big_data/transcriptomicsTorellidata/4.trimmed_reads/merged_by_species/AAL-trimmed-merged-2.fq | samtools view -Sb - | samtools sort -no - - > $TOR_PATH/big_data/transcriptomicsTorellidata/7.transcriptomeQualityCheck/B.forOneTranscriptomePerSpecies/AAL-bowtie2.nameSorted.bam

# 2.3. count 
$TRINITY_HOME/util/SAM_nameSorted_to_uniq_count_stats.pl $TOR_PATH/big_data/transcriptomicsTorellidata/7.transcriptomeQualityCheck/B.forOneTranscriptomePerSpecies/AAL-bowtie2.nameSorted.bam > AAL-uniq_count_stats

## Done for the 3 species

################ TO BE CONTINUED FROM HERE, STOPPED ON 16th OF APRIL 2021

* AAL-merged
--> longuest isoform
Stats for aligned rna-seq fragments (note, not counting those frags where neither left/right read aligned)
35944768 aligned fragments; of these:
Overall,  97.07% of aligned fragments aligned as proper pairs

* BVK-merged
Stats for aligned rna-seq fragments (note, not counting those frags where neither left/right read aligned)
46389394 aligned fragments; of these:
Overall,  97.66% of aligned fragments aligned as proper pairs

* FMN-merged
Stats for aligned rna-seq fragments (note, not counting those frags where neither left/right read aligned)
36793332 aligned fragments; of these:
Overall,  97.64% of aligned fragments aligned as proper pairs

### 3. Completeness was quantified using BUSCO (Benchmarking Universal Single-Copy Orthologs) v3.0.2 
(BUSCO: assessing genome assembly and annotation completeness with single-copy orthologs. Felipe A. Simão, Robert M. Waterhouse, Panagiotis Ioannidis, Evgenia V. Kriventseva, and Evgeny M. Zdobnov Bioinformatics, published online June 9, 2015 | Abstract | Full Text PDF | doi: 10.1093/bioinformatics/btv351) with the mamalia_odb9 database. 
BUSCO was run in mode: transcriptome

* AAL. In a total of 4104 BUSCOs, 74.9% of genes were “complete”, 5.6% were “fragmented” and 19.5% were “missing”,
considering only longuest isoforms

* BVK. In a total of 4104 BUSCOs, 75.2% of genes were “complete”, 5.3% were “fragmented” and 19.5% were “missing”,
considering only longuest isoforms

* FMN. In a total of 4104 BUSCOs, 74.0% of genes were “complete”, 4.8% were “fragmented” and 21.2% were “missing”,
considering only longuest isoforms



#### AAL-ctrl
[36,400,000]  lines read 

Stats for aligned rna-seq fragments (note, not counting those frags where neither left/right read aligned)

18245837 aligned fragments; of these:
  18245837 were paired; of these:
    559867 aligned concordantly 0 times
    17685970 aligned concordantly exactly 1 time
    0 aligned concordantly >1 times
    ----
    559867 pairs aligned concordantly 0 times; of these:
    488360 aligned as improper pairs
    71507 pairs had only one fragment end align to one or more contigs; of these:
       53957 fragments had only the left /1 read aligned; of these:
            53957 left reads mapped uniquely
            0 left reads mapped >1 times
       17550 fragments had only the right /2 read aligned; of these:
            17550 right reads mapped uniquely
            0 right reads mapped >1 times
Overall,  96.93% of aligned fragments aligned as proper pairs

## 10-3 Assesing assembly quality step 3: quantifying completeness

SOFTWARES USED :

* BUSCO v3.0.2 (BUSCO: assessing genome assembly and annotation completeness with single-copy orthologs.
Felipe A. Simão, Robert M. Waterhouse, Panagiotis Ioannidis, Evgenia V. Kriventseva, and Evgeny M. Zdobnov
Bioinformatics, published online June 9, 2015 | Abstract | Full Text PDF | doi: 10.1093/bioinformatics/btv351)

* NCBI blast v

* hmmer v.

BUSCO Benchmarking Universal Single-Copy Orthologs

Another metric of assembly quality is evaluating the extent to which it recovers single copy orthologs that are present across higher taxonomic groupings. 
While in the absence of knowing which transcripts are truly expressed in a sample it is difficult to determine an absolute expectation for recovery of these orthologs, 
clearly, high numbers of such genes classified as missing in an assembly should be considered a potential red flag. Furthermore, asssemblies based upon
the same read data can be evaluated with respect to the numbers of genes that are complete, fragmented, or missing from the assembly.

BUSCOs are selected from OrthoDB orthologous groups at major species radiations requiring orthologues to be present as single-copy genes in the vast majority (>90%) of species.
(I used mammalia)

alice@harriet:~/Dokumente/Torelli_transcriptomics/assemblyTrinityAfterTrimming/AAL-ctrl$ python /usr/local/bin/BUSCO/busco/scripts/run_BUSCO.py -i AAL-ctrl_Trinity.fasta -o AAL-ctrl_busco -l /usr/local/bin/BUSCO/busco/dataset_lineage/mamalia_odb9/ -m tran

# BUSCO version is: 3.0.2 
# The lineage dataset is: mammalia_odb9 (Creation date: 2016-02-13, number of species: 50, number of BUSCOs: 4104)
# To reproduce this run: python /usr/local/bin/BUSCO/busco/scripts/run_BUSCO.py -i AAL-ctrl_Trinity.fasta -o AAL-ctrl_busco -l /usr/local/bin/BUSCO/busco/dataset_lineage/mamalia_odb9/ -m transcriptome -c 1
#
# Summarized benchmarking in BUSCO notation for file AAL-ctrl_Trinity.fasta
# BUSCO was run in mode: transcriptome

        C:70.6%[S:12.9%,D:57.7%],F:7.2%,M:22.2%,n:4104

        2898    Complete BUSCOs (C)
        530     Complete and single-copy BUSCOs (S)
        2368    Complete and duplicated BUSCOs (D)
        297     Fragmented BUSCOs (F)
        909     Missing BUSCOs (M)
        4104    Total BUSCO groups searched

C:74.9%[S:36.9%,D:38.0%],F:5.6%,M:19.5%,n:4104

  3073  Complete BUSCOs (C)
  1514  Complete and single-copy BUSCOs (S)
  1559  Complete and duplicated BUSCOs (D)
  229   Fragmented BUSCOs (F)
  802   Missing BUSCOs (M)
  4104  Total BUSCO groups searched


Plot!
alice@harriet:~/Dokumente/Torelli_transcriptomics/assemblyTrinityAfterTrimming/AAL-ctrl$ python /usr/local/bin/BUSCO/busco/scripts/generate_plot.py -wd run_AAL-ctrl_busco/



############################################
# PART III. Look for Irg mus orthologs/paralogs by tblastn 

# Extraction of IRGb2b1 for the 3 transcriptomes, by homologies with different sequences known
blastn -query ~/Dokumente/Torelli_transcriptomics/blastTandems/Tandems.fa 
       -db ~/Dokumente/Torelli_transcriptomics/assemblyOneTranscriptomPerSpecies/AAL-merged_trinity/AAL-merged-longuest-isoform-transcriptome.fasta 
       -outfmt 11 
       -evalue 1e-50 
       -out AAL-mlit.vs.tandems.blastn

And select the uniq transcripts
--> too fragmented....

############################################
# PART IV. RNAseq analysis

# Strategy
# * map reads back to assembly and count them
# * DE (differential expression) analysis between CTRL and IFNg


# Building count matrix for the 6 samples:
$TRINITY_HOME/util/align_and_estimate_abundance.pl --transcripts ../AAL-Trinity-transcriptome.fasta --seqType fq --left ../../../rawData/filtratedRNAreads/AAL-ctrl_1.cor.fq --right ../../../rawData/filtratedRNAreads/AAL-ctrl_2.cor.fq --est_method RSEM --output_dir AAL-ctrl --aln_method bowtie2 --prep_reference

#################
## EXTRA STEP: full RNAseq analysis (to complete)

# GO, KEGG, "name", whatever, Rbiomat (uniref/uniprot) after DEseq

# For orthofinder: extract longuest isoform per trinity gene using get_longest_isoform_seq_per_trinity_gene.pl script from Trinity.
# source https://github.com/davidemms/OrthoFinder/issues/38 FOR ORTHOFINDER

## Detection of orthologous genes between species
# Transcriptomes were translated in predicted proteins files using TransDecoder v5.3.0 (TransDecoder.LongOrfs and TransDecoder.Predict function)
# Reference proteome of Mus musculus (mouse) (Strain: C57BL/6J) was used for further annotation 
# (Proteome ID: P000000589, Organism ID: 10090, https://www.uniprot.org/proteomes/, 20,210 predicted protein-coding genes. 
# From GRCm38.p6 (Genome Reference Consortium Mouse Reference 38), INSDC Assembly GCA_000001635.8, Jan 2012)

#We used Orthofinder (Emms D.M. & Kelly S. (2015), Genome Biology 16:157) on our four translated transcriptomes and our Mus musculus reference.
#An orthogroup is defined as the group of genes descended from a single gene in the last common ancestor of a group of species.
#We want to extract orthologous genes, therefore we focus on single-copy orthologous genes to exclude potential paralogues.
#OrthoFinder assigned 67557 genes (73.7% of total) to 13728 orthogroups. 
#Half of the assigned genes were in orthogroups of 5 genes or more (G50 = 5).
#There were 8725 orthogroups wit6h all species present and 1367 of these consisted entirely of single-copy genes.
#We further focus on these 1367 genes, orthologues amongst all considered rodents.

#* then you can retrieve nucleotide sequences of OrthologsIDS.txt from species`s nucleotide files.
#grep -Fwf OrthoFinder/Results_Dec23/Orthogroups/Orthogroups_SingleCopyOrthologues.txt OrthoFinder/Results_Dec23/Orthogroups/Orthogroups.txt > OrthoFinder/Results_Dec23/OrthologsIDS.txt
#and grep in cds files created by transdecoder, corresponding to predicted proteins 
#(seqtk subseq ApodemusIds.txt cds_files/Apodemus_agrarius.cds)

#--> FinalTranscriptomes_onlyOrthologous

#Finally, we uniformise the fasta headers to Mus orhologous ones

#https://www.uniprot.org/uploadlists/ get gene ensemble names from uniprotKB names

#table chrAliases 
#comma delimited text file that includes two columns.
#The first column gives the chromosome names used in the annotation and the second column gives the chromosome names used by reads

#awk -v OFS="\t" '$1=$1' orthologous.list | cut -f2,5 > temp
#sed 's/\t/,/g' temp > Apodemus_agrarius_chrAliases.csv

#If you have multiple RNA-Seq data sets that you want to compare (eg. different tissues sampled from a single organism), be sure to generate a single Trinity assembly and to then run the abundance estimation separately for each of your samples.

#then, align with subread (v1.6.3) the 6 sets of reads to the 3 transcriptomes

#- 1. build index
#subread-1.6.3-Linux-x86_64/bin/subread-buildindex -o myindexApodemus ../orthologousSearch/Apodemus_agrarius_OG.csd 

#- 2. align
#FeatureCountAnalysis$ 
#subread-1.6.3-Linux-x86_64/bin/subread-align -i myindexApodemus -r ../rawData/trimmed_reads/fixed_AAL-ctrl_1.cor_val_1.fq -R ../rawData/trimmed_reads/fixed_AAL-ctrl_2.cor_val_2.fq -t 0 -o ResultsApodemus_ctrl.bam

#- 3. Counted featureCounts
#FeatureCountAnalysis$ 
#subread-1.6.3-Linux-x86_64/bin/featureCounts -a gtf_files/Apodemus_agrarius_OG.fasta.gtf -p -t exon -o Apodemus_agrarius_ctrl_FC.txt -A ResultsApodemus_ctrl.bam 

#NB error chromosome identifier between alignment in bam and GTF
#samtools view -h -o out.sam in.bam

#*** make GTF ***
#Torelli_transcriptomics/orthologousSearch/test$ 
#blat ../FinalTranscriptomes_onlyOrthologous/Apodemus_agrarius_OG.fasta ../../FeatureCountAnalysis/GCF_000001635.26_GRCm38.p6_genomic.fna temp.psl -t=dna -q=dna -tileSize=11 -minIdentity=90 -maxIntron=1001 -out=psl

#Idea. in GTF, replace first column chromosome by the Mus proteome identifier
#Plus simple. Par la colonnu gene, et pas de chrAlias

#Idea bonus! 

#lancer un blastx vs trembl et swissprot des 3 fulls transcritpomes pour annotation. Ca va prendre un bail.

#/blastxForAnnotation$ for FILE in ../assemblyOneTranscriptomPerSpecies/*/*transcriptome.fasta ; do bash script.sh $FILE; done 

#!/bin/sh/bash

# Done before: makeblastdb -in uniprot_sprot.fasta -dbtype prot; makeblastdb -in uniprot_sprot.fasta -dbtype prot

#fbname=$(basename "$1" .txt)

#echo "$fbname"

#blastx -query $1 -db uniprot_sprot.fasta -num_threads 8 -max_target_seqs 1 -outfmt 6 -out blastx_sprot_$fbname.outfmt6

#echo "$fbname"

#blastx -query $1 -db uniprot_sprot.fasta -num_threads 8 -max_target_seqs 1 -outfmt 6 -out blastx_trembl_$fbname.outfmt6

#TO DO : look at the sequences of the found IRG

# DE analysis

#WHY featureCount/edgeR: https://www.researchgate.net/publication/328753474_A_comprehensive_RNA-Seq_pipeline_includes_meta-analysis_interactivity_and_automatic_reporting
#"featureCounts-edgeRpipelines are the ones that give better results in differential analysis, as most of the genes they identifyare in common with other methods"

#We performed a GO analysis on the subsets of transcripts that were identified by OrthoFinder as being single copy orthologous between the four species (5 samples)

#subread to align to homologous single copy
#feature count for coverage
