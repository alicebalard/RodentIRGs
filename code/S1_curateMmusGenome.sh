#!/bin/bash

cd /SAN/Alices_sandpit/Torelli_transcriptomics/big_data/02proteomes/14rodents

cp Mmus_GCF_000001635.27_GRCm39_protein.faa Mmus_GCF_000001635.27_GRCm39_CURATED_protein.faa
## replace IRGs proteins in the curated Mmus proteome by correctly curated sequences

## 1. transform the proteome in one liner for further grep -A
awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' Mmus_GCF_000001635.27_GRCm39_CURATED_protein.faa > temp
mv temp Mmus_GCF_000001635.27_GRCm39_CURATED_protein.faa # 184997 lines

## 2. check proteins sequences
irgSeqs="/SAN/Alices_sandpit/Torelli_transcriptomics/GIT/RodentIRGs/data/fasta_sequences/ReferenceMMus_codingIrgs_CURATED.fasta"
curatGen="/SAN/Alices_sandpit/Torelli_transcriptomics/big_data/02proteomes/14rodents/Mmus_GCF_000001635.27_GRCm39_CURATED_protein.faa"

## Irg that just need a change of name:
# "NP_001028939.1" irga2 ; "NP_001030031.2" irga3 ; "NP_001094945.1" irga4 ; "NP_068564.4" irga6 ; "NP_001342686.1" irgm1 ; "NP_001039005.1" irgb2b1
# "NP_001258606.1" irgd ; "NP_001128587.1" irgb10 ; "NP_061208.3" irgm3 ; "NP_062313.3" irgm2 ; "NP_001344963.1" irgc ; "NP_694774.3" irgq
# "NP_001138636.1" irgb6
# # >Mmus_Irgb6* NP_035709.3 check that it's removed

printf "NP_001028939.1\nNP_001030031.2\nNP_001094945.1\nNP_068564.4\nNP_001342686.1\nNP_001039005.1\nNP_001258606.1\nNP_001128587.1\nNP_061208.3\nNP_062313.3\nNP_001344963.1\nNP_694774.3\nNP_001138636.1" > temp_irgList1

cat temp_irgList1 | while read line
do
    grep -A 1 $line $irgSeqs > tempIrgFasta1 # store the full IRG's new fasta
    seq1=$(grep -A 1 $line $irgSeqs | sed -n '2p') # only the irg sequence
    name2=$(grep -A 1 $line $curatGen | sed -n '1p') # irg name in original genome
    seq2=$(grep -A 1 $line $curatGen | sed -n '2p') # irg seq in original genome
    ## compare sequences to make sure they are equal
    if [ "$seq1" = "$seq2" ]; then
	echo "Strings are equal."
	## 3. remove multiple transcripts from curated genome 
	grep -B 1 $seq2 $curatGen | grep -v '^--' > temp
	awk 'NR==FNR{a[$0];next} !($0 in a)' temp $curatGen > tempGenome
	## 4. add new curated IRG sequence
	cat tempIrgFasta1 tempGenome > Mmus_GCF_000001635.27_GRCm39_CURATED_protein.faa
    else
	echo "Strings are not equal."
	## 4. only add new curated IRG sequence
	cat tempIrgFasta1 Mmus_GCF_000001635.27_GRCm39_CURATED_protein.faa > tempGenome
	mv tempGenome Mmus_GCF_000001635.27_GRCm39_CURATED_protein.faa
    fi
done

## new curations (done one by one):

#### "CBY65983.1" irga7
grep -A 1 "Mmus_Irga7" $irgSeqs > tempIrgFasta_toadd # store the full IRG's new fasta
## add new curated IRG sequence
cat tempIrgFasta_toadd tempGenome > Mmus_GCF_000001635.27_GRCm39_CURATED_protein.faa

#### b9b8
grep -A 1 "NEW_CURATED_Irgb9b8" $irgSeqs > tempIrgFasta_toadd # store the full IRG's new fasta
grep -A 1 "NP_001019401.2" Mmus_GCF_000001635.27_GRCm39_CURATED_protein.faa > tempIrgFasta_torm # the wrong irg to rm
seq2=$(cat tempIrgFasta_torm | sed -n '2p') # irg seq in original genome
## remove multiple transcripts from curated genome 
grep -B 1 $seq2 $curatGen | grep -v '^--' > temp
awk 'NR==FNR{a[$0];next} !($0 in a)' temp $curatGen > tempGenome
## add new curated IRG sequence
cat tempIrgFasta_toadd tempGenome > Mmus_GCF_000001635.27_GRCm39_CURATED_protein.faa

#### b5*b3
grep -A 1 "NEW_CURATED_Irgb5\*b3" $irgSeqs > tempIrgFasta_toadd # store the full IRG's new fasta
grep -A 1 "NP_001108151.1" Mmus_GCF_000001635.27_GRCm39_CURATED_protein.faa > tempIrgFasta_torm # the wrong irg to rm
seq2=$(cat tempIrgFasta_torm | sed -n '2p') # irg seq in original genome
## remove multiple transcripts from curated genome 
grep -B 1 $seq2 $curatGen | grep -v '^--' > temp
awk 'NR==FNR{a[$0];next} !($0 in a)' temp $curatGen > tempGenome
## add new curated IRG sequence
cat tempIrgFasta_toadd tempGenome > Mmus_GCF_000001635.27_GRCm39_CURATED_protein.faa

#### b5b4
grep -A 1 "NEW_CURATED_Irgb5b4" $irgSeqs > tempIrgFasta_toadd # store the full IRG's new fasta
grep -A 1 "NP_775610.1" Mmus_GCF_000001635.27_GRCm39_CURATED_protein.faa > tempIrgFasta_torm # the wrong irg to rm
seq2=$(cat tempIrgFasta_torm | sed -n '2p') # irg seq in original genome
## remove multiple transcripts from curated genome 
grep -B 1 $seq2 $curatGen | grep -v '^--' > temp
awk 'NR==FNR{a[$0];next} !($0 in a)' temp $curatGen > tempGenome
## add new curated IRG sequence
cat tempIrgFasta_toadd tempGenome > Mmus_GCF_000001635.27_GRCm39_CURATED_protein.faa

## remove all temporary files
rm temp*
