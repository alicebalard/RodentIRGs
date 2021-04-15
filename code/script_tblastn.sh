# A. Balard
# 13 April 2021
## all done on Harriet server from Heitlinger team
# path for data (on Harriet): 
TOR_PATH=/SAN/Alices_sandpit/Torelli_transcriptomics

## A. Identify nucleotide and protein sequences from Bekpen to the newest Mus musculus genome
tblastn -db $TOR_PATH/big_data/genomes/musMusculus/GCF_000001635.27_GRCm39_genomic.fna -query $TOR_PATH/GIT/RodentIRGs/data/REFERENCE_Protein_sequences_of_IRG_mouse_GRCm39_and_Bekpen2005 -num_threads 16 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sseq" -out $TOR_PATH/GIT/RodentIRGs/data/blast_results/IRGmus2mus/tblastnIRGBekpen_vs_MmusGRCm39.outfmt6 -max_target_seqs 3 -max_hsps 2

## B. Search for IRGs: TBLASTN search translated nucleotide databases using a protein query in genomes for our 11 rodent species
for F in $TOR_PATH/big_data/genomes/*.fna; do makeblastdb -in $F -dbtype nucl; done
for F in $TOR_PATH/big_data/genomes/*.fna; do tblastn -query $TOR_PATH/GIT/RodentIRGs/data/REFERENCE_Protein_sequences_of_IRG_mouse_GRCm39_and_Bekpen2005.fasta -db $F -num_threads 16 -evalue 1e-100 -max_hsps 10 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sseq" -out $TOR_PATH/GIT/RodentIRGs/data/blast_results/IRGmus2otherRodents/tBlastnIRGmusGRCm39_vs_$(basename -- "$F").outfmt6; done

## Extra tools:

## to find specific bits in a fasta genome:
# samtools faidx big_data/genomes/musMusculus/GCF_000001635.27_GRCm39_genomic.fna NC_000084.7:60574825-60576058
## to get reverse complement:
# samtools faidx big_data/genomes/musMusculus/GCF_000001635.27_GRCm39_genomic.fna NC_000077.7:48909897-48911170 | seqtk seq -r

## transform in fasta
# for myfile in *outfmt6; do awk 'BEGIN { OFS = "\n" } { print ">"substr(FILENAME,1,4)"_"$1" "$2":"$9"-"$10" evalue:"$11" mismatches:"$5, $13 }' $myfile > $myfile.fasta; done

## concatenate
# for i in *.fasta; do cat $i >> concatenatedtblastn.fasta; done
