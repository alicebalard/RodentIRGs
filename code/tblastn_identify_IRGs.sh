## IRG identification by tblasn
## A. Balard, 14 April 2021

# path for data (on Harriet): /SAN/Alices_sandpit/Torelli_transcriptomics/

## A. Identify nucleotide and protein sequences from Bekpen to the newest Mus musculus genome
tblastn -db PATH/big_data/genomes/musMusculus/GCF_000001635.27_GRCm39_genomic.fna -query PATH/data/Protein_sequences_of_IRG_mouse_Bekpen_2005.fasta -num_threads 16 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sseq" -out ../data/tblastnIRGBekpen_vs_MmusGRCm39.outfmt6 -max_target_seqs 3 -max_hsps 2
blastn -db PATH/big_data/genomes/musMusculus/GCF_000001635.27_GRCm39_genomic.fna -query PATH/data/Nucleotide_sequences_of_IRG_mouse_Bekpen_2005.fasta -num_threads 16 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sseq" -out ../data/blastnIRGBekpen_vs_MmusGRCm39.outfmt6 -max_target_seqs 3 -max_hsps 2

## to find specific bits in a fasta genome (to complete the sequences):
# samtools faidx big_data/genomes/musMusculus/GCF_000001635.27_GRCm39_genomic.fna NC_000084.7:60574825-60576058
## to get reverse complement:
# samtools faidx big_data/genomes/musMusculus/GCF_000001635.27_GRCm39_genomic.fna NC_000077.7:48909897-48911170 | seqtk seq -r

## B. Search for IRGs: TBLASTN search translated nucleotide databases using a protein query in genomes for our 11 rodent species
for F in PATH/big_data/genomes/*.fna; do makeblastdb -in $F -dbtype nucl; done
for F in PATH/big_data/genomes/*.fna; do tblastn -query PATH/GIT/RodentIRGs/figures/Supplementary/S2_Protein_sequences_of_IRG_mouse_GRCm39.fasta -db $F -num_threads 16 -evalue 1e-100 -max_hsps 10 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sseq" -out PATH/GIT/RodentIRGs/data/blast_results/IRGmus2otherRodents/tBlastnIRGmusGRCm39_vs_$(basename -- "$F").outfmt6; done

## transform in fasta
# for myfile in *outfmt6; do awk 'BEGIN { OFS = "\n" } { print ">"substr(FILENAME,1,4)"_"$1" "$2":"$9"-"$10" evalue:"$11" mismatches:"$5, $13 }' $myfile > $myfile.fasta; done

## concatenate
# for i in *.fasta; do cat $i >> concatenatedtblastn.fasta; done

## C. Run orthofinder
### The files from Ensembl will contain many transcripts per gene. If we ran OrthoFinder on these raw files it would take 10x longer than necessary and could lower the accuracy. We’ll use a script provided with OrthoFinder to extract just the longest transcript variant per gene and run OrthoFinder on these files:
# for f in *faa ; do python /localstorage/alice/anaconda3/pkgs/orthofinder-2.5.1-0/bin/primary_transcript.py $f ; done
### Run Orthofinder
# orthofinder -f primary_transcripts/

# Extract gene trees for IRGs
## from GIT/RodentIRGs$ 
cat figures/Table1_IRGs_naming.csv | cut -d ',' -f6 > tempIRGprot.txt;
cat tempIRGprot.txt | sed '/^$/d' > temp2;
cat temp2 | sed '/?/d' > tempIRGprot.txt;
sed 's/\;/\n/g' tempIRGprot.txt > temp;
sed 's/ //g' temp > tempIRGprot.txt;
rm temp; rm temp2

# Extract all orthologues: 
TOR_PATH=/SAN/Alices_sandpit/Torelli_transcriptomics/;
head -n 1 $TOR_PATH/big_data/proteomes/Results_OrthoFinder/Phylogenetic_Hierarchical_Orthogroups/N0.tsv > $TOR_PATH/GIT/RodentIRGs/rawdata/orthoFinder_sub/N0_HOG_IRG.tsv;
rg -f $TOR_PATH/GIT/RodentIRGs/tempIRGprot.txt $TOR_PATH/big_data/proteomes/Results_OrthoFinder/Phylogenetic_Hierarchical_Orthogroups/N0.tsv >> $TOR_PATH/GIT/RodentIRGs/rawdata/orthoFinder_sub/N0_HOG_IRG.tsv


