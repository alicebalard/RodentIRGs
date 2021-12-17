# A. Balard
# October 2021
# Analyses done on Harriet server from Heitlinger team

# path for data (on Harriet):
TOR_PATH=/SAN/Alices_sandpit/Torelli_transcriptomics

## A. Identify nucleotide and protein sequences from Bekpen to the newest Mus musculus genome
makeblastdb -in $TOR_PATH/big_data/genomes/Mmus_GCF_000001635.27_GRCm39_genomic.fna -dbtype "nucl"
tblastn -db $TOR_PATH/big_data/genomes/Mmus_GCF_000001635.27_GRCm39_genomic.fna -query $TOR_PATH/GIT/RodentIRGs/data/REFERENCE_Protein_sequences_of_IRG_mouse_GRCm39_and_Bekpen2005.fasta -num_threads 16 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sseq" -out $TOR_PATH/GIT/RodentIRGs/data/blast_results/IRGmus2mus/tblastnIRGBekpen_vs_MmusGRCm39.outfmt6 -max_target_seqs 3 -max_hsps 2

#################
## update december 2021: replace proteome Mmus with curated CORRECT IRG sequences
### Sequences in: RodentIRGs/data/fasta_sequences/ReferenceMMus_codingIrgs_CURATED.fasta 

bash /SAN/Alices_sandpit/Torelli_transcriptomics/GIT/RodentIRGs/code/S1_curateMmusGenome.sh

################
## B. OrthoFinder identification
## source: https://github.com/davidemms/OrthoFinder
### The files from Ensembl will contain many transcripts per gene. If we ran OrthoFinder on these raw files it would take 10x longer than necessary and could lower the accuracy. We’ll use a script provided with OrthoFinder to extract just the longest transcript variant per gene and run OrthoFinder on these files:
cd /SAN/Alices_sandpit/Torelli_transcriptomics/big_data/02proteomes/14rodents/
for f in *faa ; do python /localstorage/alice/anaconda3/pkgs/orthofinder-2.5.1-0/bin/primary_transcript.py $f ; done

### Run Orthofinder
orthofinder -f primary_transcripts/

# Extract IRGs Mmus protein codes 
grep ">" /SAN/Alices_sandpit/Torelli_transcriptomics/GIT/RodentIRGs/data/fasta_sequences/ReferenceMMus_codingIrgs_CURATED.fasta | cut -d ' ' -f1 | sed 's/>//' > /SAN/Alices_sandpit/Torelli_transcriptomics/GIT/RodentIRGs/data/orthofinder_results/tempIRGprot.txt

## I'm here




# Extract all orthologues:

## headers
head -n 1 $TOR_PATH/big_data/02proteomes/14rodents/primary_transcripts/OrthoFinder/Results_Sep16_1/Phylogenetic_Hierarchical_Orthogroups/N0.tsv > $TOR_PATH/GIT/RodentIRGs/data/orthofinder_results/N0_HOG_IRG.tsv;
## OF: N0.tsv is a tab separated text file. Each row contains the genes belonging to a single orthogroup. The genes from each orthogroup are organized into columns, one per species. Additional columns give the HOG (Hierarchical Orthogroup) ID and the node in the gene tree from which the HOG was determined (note, this can be above the root of the clade containing the genes). This file effectively replaces the orthogroups in Orthogroups/Orthogroups.tsv from Markov clustering using MCL.
rg -f $TOR_PATH/GIT/RodentIRGs/data/orthofinder_results/tempIRGprot.txt $TOR_PATH/big_data/02proteomes/14rodents/primary_transcripts/OrthoFinder/Results_Sep16_1/Phylogenetic_Hierarchical_Orthogroups/N0.tsv >> $TOR_PATH/GIT/RodentIRGs/data/orthofinder_results/N0_HOG_IRG.tsv
rm $TOR_PATH/GIT/RodentIRGs/data/orthofinder_results/tempIRGprot.txt

# We match the OrthoFinder results with the known IRG names in Mmus and extract the fasta sequences
## See mainScript.R

## C. Reciprocal Best Blast Hit identification:

# Extract IRGs Mmus protein codes 
cat $TOR_PATH/GIT/RodentIRGs/data/Table1_IRGs_naming.csv | cut -d ',' -f1,6 | sed 's/\;.*//' | sed '/?/d' | awk -F, '$2 != ""' > $TOR_PATH/GIT/RodentIRGs/data/GRCm39IRGprotseq.txt

# Extract IRGs Mmus protein sequences
######
# -->done in R mainScript.R
# $TOR_PATH/GIT/RodentIRGs/data/GRCm39IRGprotseq.fasta
######

# make blast db from our additional rodent genomes
for F in $TOR_PATH/big_data/01genomes/3extraRodents/*.fna; do makeblastdb -in $F -dbtype nucl; done

# make blast db from the Mmus IRG (original query)
makeblastdb -in $TOR_PATH/GIT/RodentIRGs/data/fasta_sequences/Mmus_IRG_protseq.fasta -dbtype prot

## Species one by one
#Input original query file
mouseIRGquery1=$TOR_PATH/GIT/RodentIRGs/data/fasta_sequences/Mmus_IRG_protseq.fasta
#Path to output results
outputPath=$TOR_PATH/GIT/RodentIRGs/data/blast_results/RBBS_Irgmus2set2

# first blast search: tblastn mus IRG protein sequences vs 3 extra rodent genomes (4 as one has 2)
for F in $TOR_PATH/big_data/01genomes/3extraRodents/*.fna; do tblastn -query $mouseIRGquery1 -db $F -num_threads 60 -evalue 1e-100 -max_hsps 1 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sseq" -out $outputPath/tBlastnIRGmusGRCm39_vs_$(basename -- "$F").outfmt6; done

# # convert blast results into fasta 
# for myfile in $outputPath/*outfmt6; do
#     newname=$(echo "${myfile##*/}" | sed -n -e 's/^.*tBlastnIRGmusGRCm39_vs_//p' |  awk '{print substr($0, 1, 4);exit}') ;
#     awk -v var="$newname" 'BEGIN { OFS = "\n" } { print ">" var "_"$1" "$2":"$9"-"$10" evalue:"$11" mismatches:"$5, $13 }' $myfile > $myfile.fasta;
# done

# reciprocal BLAST search: BLASTX search protein databases using a translated nucleotide query
for F in $TOR_PATH/big_data/01genomes/3extraRodents/*.fna; do blastx -db $mouseIRGquery1 -query $F -num_threads 60 -evalue 1e-100 -max_hsps 1 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sseq" -out $outputPath/tBlastn$(basename -- "$F")_vs_IRGmusGRCm39.outfmt6; done

# # convert blast results into fasta 
# for myfile in $outputPath/*vs_IRGmusGRCm39.outfmt6; do
#     newname=$(echo "${myfile##*/}" | sed -n -e 's/^.*tBlastn//p' |  awk '{print substr($0, 1, 4);exit}') ;
#     awk -v var="$newname" 'BEGIN { OFS = "\n" } { print ">" var "_"$1" "$2":"$9"-"$10" evalue:"$11" mismatches:"$5, $13 }' $myfile > $myfile.fasta;
# done

# Load the blast results into R (see mainScript.R)
for file in *; do cat $file | sed -e 's/\t/ /g' > $file.2; done
## This is done for easier handling in R 
