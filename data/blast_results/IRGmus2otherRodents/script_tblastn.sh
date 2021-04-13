# A. Balard
# 13 April 2021

## tBlastn mus IRGs in each rodent genomes (N=11 species)
for F in /SAN/Alices_sandpit/Torelli_transcriptomics/big_data/genomes/*.fna; 
do tblastn -query /SAN/Alices_sandpit/Torelli_transcriptomics/GIT/RodentIRGs/figures/Supplementary/S2_Protein_sequences_of_IRG_mouse_GRCm39.fasta -db $F -num_threads 16 -evalue 1e-100 -max_hsps 10 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sseq" -out /SAN/Alices_sandpit/Torelli_transcriptomics/GIT/RodentIRGs/data/blast_results/IRGmus2otherRodents/tBlastnIRGmusGRCm39_vs_$(basename -- "$F").outfmt6; done

## transform blast hits in fasta within each species
for myfile in *outfmt6; do awk 'BEGIN { OFS = "\n" } { print ">"substr(FILENAME,24,27)"_"$1" "$2":"$9"-"$10" evalue:"$11" mismatches:"$5, $13 }' $myfile > $myfile.fasta; done

## remove duplicates

## blast in the other sense: putativeIRG searched in mouse IRG list
