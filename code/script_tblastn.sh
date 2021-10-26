# A. Balard
# 13 April 2021
## all done on Harriet server from Heitlinger team
# path for data (on Harriet): 
TOR_PATH=/SAN/Alices_sandpit/Torelli_transcriptomics

## A. Identify nucleotide and protein sequences from Bekpen to the newest Mus musculus genome
makeblastdb -in $TOR_PATH/big_data/genomes/Mmus_GCF_000001635.27_GRCm39_genomic.fna -dbtype "nucl"
tblastn -db $TOR_PATH/big_data/genomes/Mmus_GCF_000001635.27_GRCm39_genomic.fna -query $TOR_PATH/GIT/RodentIRGs/data/REFERENCE_Protein_sequences_of_IRG_mouse_GRCm39_and_Bekpen2005.fasta -num_threads 16 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sseq" -out $TOR_PATH/GIT/RodentIRGs/data/blast_results/IRGmus2mus/tblastnIRGBekpen_vs_MmusGRCm39.outfmt6 -max_target_seqs 3 -max_hsps 2

## B. Search for IRGs: TBLASTN search translated nucleotide databases using a protein query in genomes for our 11 rodent species
for F in $TOR_PATH/big_data/genomes/*.fna; do makeblastdb -in $F -dbtype nucl; done
for F in $TOR_PATH/big_data/genomes/*.fna; do tblastn -query $TOR_PATH/GIT/RodentIRGs/data/REFERENCE_Protein_sequences_of_IRG_mouse_GRCm39_and_Bekpen2005.fasta -db $F -num_threads 16 -evalue 1e-100 -max_hsps 10 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sseq" -out $TOR_PATH/GIT/RodentIRGs/data/blast_results/IRGmus2otherRodents/tBlastnIRGmusGRCm39_vs_$(basename -- "$F").outfmt6; done

## C. search in transcriptomes
T1=/SAN/Alices_sandpit/Torelli_transcriptomics/big_data/transcriptomicsTorellidata/6.TrinityAssemblies/B.oneTranscriptomePerSpecies/AAL-merged_trinity/AAL-Trinity-transcriptome.fasta
T2=/SAN/Alices_sandpit/Torelli_transcriptomics/big_data/transcriptomicsTorellidata/6.TrinityAssemblies/B.oneTranscriptomePerSpecies/BVK-merged_trinity/BVK-Trinity-transcriptome.fasta
T3=/SAN/Alices_sandpit/Torelli_transcriptomics/big_data/transcriptomicsTorellidata/6.TrinityAssemblies/B.oneTranscriptomePerSpecies/FMN-merged_trinity/FMN-Trinity-transcriptome.fasta

for f in $T1 $T2 $T3; do makeblastdb -in $f -dbtype nucl; done
for f in $T1 $T2 $T3; do tblastn -query $TOR_PATH/GIT/RodentIRGs/data/REFERENCE_Protein_sequences_of_IRG_mouse_GRCm39_and_Bekpen2005.fasta -db $f -num_threads 16 -evalue 1e-100 -max_hsps 10 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sseq" -out $TOR_PATH/GIT/RodentIRGs/data/blast_results/IRGmus2transcriptomes/tBlastnIRGmusGRCm39_vs_$(basename -- "$f").outfmt6; done

## In GG transcriptome
T1=/SAN/Alices_sandpit/Torelli_transcriptomics/big_data/transcriptomicsTorellidata/6.TrinityAssemblies/C.guidedTrinity/AAL-trinity_out_dir/AAL-Trinity-GG.fasta # rename!
T2=/SAN/Alices_sandpit/Torelli_transcriptomics/big_data/transcriptomicsTorellidata/6.TrinityAssemblies/C.guidedTrinity/BVK-trinity_out_dir/BVK-Trinity-GG.fasta # rename!
T3=/SAN/Alices_sandpit/Torelli_transcriptomics/big_data/transcriptomicsTorellidata/6.TrinityAssemblies/C.guidedTrinity/FMN-trinity_out_dir/FMN-Trinity-GG.fasta # rename!

for f in $T1 $T2 $T3; do makeblastdb -in $f -dbtype nucl; done
for f in $T1 $T2 $T3; do tblastn -query $TOR_PATH/GIT/RodentIRGs/data/REFERENCE_Protein_sequences_of_IRG_mouse_GRCm39_and_Bekpen2005.fasta -db $f -num_threads 16 -evalue 1e-100 -max_hsps 10 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sseq" -out $TOR_PATH/GIT/RodentIRGs/data/blast_results/IRGmus2transcriptomes/tBlastnIRGmusGRCm39_vs_$(basename -- "$f").outfmt6; done

## Extra tools:

## to find specific bits in a fasta genome:
# samtools faidx big_data/genomes/musMusculus/GCF_000001635.27_GRCm39_genomic.fna NC_000084.7:60574825-60576058
## to get reverse complement:
# samtools faidx big_data/genomes/musMusculus/GCF_000001635.27_GRCm39_genomic.fna NC_000077.7:48909897-48911170 | seqtk seq -r

## transform in fasta
# for myfile in *outfmt6; do awk 'BEGIN { OFS = "\n" } { print ">"substr(FILENAME,1,4)"_"$1" "$2":"$9"-"$10" evalue:"$11" mismatches:"$5, $13 }' $myfile > $myfile.fasta; done

## concatenate
# for i in *.fasta; do cat $i >> concatenatedtblastn.fasta; done
