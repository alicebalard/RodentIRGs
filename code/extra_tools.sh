# A. Balard

## to find specific bits in a fasta genome:
# samtools faidx big_data/genomes/musMusculus/GCF_000001635.27_GRCm39_genomic.fna NC_000084.7:60574825-60576058
## to get reverse complement:
# samtools faidx big_data/genomes/musMusculus/GCF_000001635.27_GRCm39_genomic.fna NC_000077.7:48909897-48911170 | seqtk seq -r

## transform in fasta
# for myfile in *outfmt6; do awk 'BEGIN { OFS = "\n" } { print ">"substr(FILENAME,1,4)"_"$1" "$2":"$9"-"$10" evalue:"$11" mismatches:"$5, $13 }' $myfile > $myfile.fasta; done

## concatenate
# for i in *.fasta; do cat $i >> concatenatedtblastn.fasta; done
