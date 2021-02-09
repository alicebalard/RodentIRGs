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

## C. Run orthofinder
### The files from Ensembl will contain many transcripts per gene. If we ran OrthoFinder on these raw files it would take 10x longer than necessary and could lower the accuracy. We’ll use a script provided with OrthoFinder to extract just the longest transcript variant per gene and run OrthoFinder on these files:
# for f in *faa ; do python /localstorage/alice/anaconda3/pkgs/orthofinder-2.5.1-0/bin/primary_transcript.py $f ; done
### Run Orthofinder
# orthofinder -f primary_transcripts/

## extract orthogroups
# grep -E "NP_001028939.1|NP_001030031.2|NP_001094945.1|NP_068564.4|NP_001342686.1|NP_062313.3|NP_061208.3|NP_001039005.1|NP_001019401.2|NP_001108151.1|NP_035709.3|NP_775610.1|NP_001138636.1|NP_001258606.1|NP_001128587.1|NP_001344963.1|NP_694774.3" primary_transcripts/OrthoFinder/Results_Feb06/Orthogroups/Orthogroups.tsv > Mmus_GRCm39_IRG_prot_names_corresp_orthofinder.csv
## -> gitted file
