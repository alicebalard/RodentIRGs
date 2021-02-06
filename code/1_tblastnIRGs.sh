## A. Identify nucleotide and protein sequences from Bekpen to the newest Mus musculus genome
# tblastn -db ../../big_data/genomes/musMusculus/GCF_000001635.27_GRCm39_genomic.fna -query ../data/Protein_sequences_of_IRG_mouse_Bekpen_2005.fasta -num_threads 16 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sseq" -out ../data/tblastnIRGBekpen_vs_MmusGRCm39.outfmt6 -max_target_seqs 3 -max_hsps 2
# blastn -db ../../big_data/genomes/musMusculus/GCF_000001635.27_GRCm39_genomic.fna -query ../data/Nucleotide_sequences_of_IRG_mouse_Bekpen_2005.fasta -num_threads 16 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sseq" -out ../data/blastnIRGBekpen_vs_MmusGRCm39.outfmt6 -max_target_seqs 3 -max_hsps 2

## to find specific bits in a fasta genome (to complete the sequences):
# samtools faidx big_data/genomes/musMusculus/GCF_000001635.27_GRCm39_genomic.fna NC_000084.7:60574825-60576058
## to get reverse complement:
# samtools faidx big_data/genomes/musMusculus/GCF_000001635.27_GRCm39_genomic.fna NC_000077.7:48909897-48911170 | seqtk seq -r

## B. Search for IRGs: TBLASTN search translated nucleotide databases using a protein query in genomes for our 9 rodent species 
# for F in /SAN/Alices_sandpit/Torelli_transcriptomics/big_data/genomes/*.fna; do makeblastdb -in $F -dbtype nucl; done
# for F in /SAN/Alices_sandpit/Torelli_transcriptomics/big_data/genomes/*.fna; do tblastn -query ../data/fasta/Protein_sequences_of_IRG_mouse_GRCm39.fasta -db $F -num_threads 16 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sseq" -out ../data/blast_results/tblstnIRGmusGRCm39_vs_${F%.*}.outfmt6

## move results in good folder
# for i in *outfmt6; do mv -- "$i" "${i/%.outfmt6/_queryIRG_GRCm39.outfmt6}"; done
# moved in: RodentIRGs/data/blast_results/IRGmus2otherRodents

## transform in fasta
# for myfile in *outfmt6; do awk 'BEGIN { OFS = "\n" } { print ">"substr(FILENAME,1,4)"_"$1" "$2":"$9"-"$10" evalue:"$11" mismatches:"$5, $13 }' $myfile > $myfile.fasta; done

## concatenate
# for i in *.fasta; do cat $i >> concatenatedtblastn.fasta; done

