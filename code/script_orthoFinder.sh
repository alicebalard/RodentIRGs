## OrthoFinder analysis script
## source: https://github.com/davidemms/OrthoFinder
## all done on Harriet server from Heitlinger team

### The files from Ensembl will contain many transcripts per gene. If we ran OrthoFinder on these raw files it would take 10x longer than necessary and could lower the accuracy. We’ll use a script provided with OrthoFinder to extract just the longest transcript variant per gene and run OrthoFinder on these files:

for f in *faa ; do python /localstorage/alice/anaconda3/pkgs/orthofinder-2.5.1-0/bin/primary_transcript.py $f ; done

### Run Orthofinder

orthofinder -f primary_transcripts/

# Extract gene trees for IRGs
TOR_PATH=/SAN/Alices_sandpit/Torelli_transcriptomics;

cat $TOR_PATH/GIT/RodentIRGs/data/Table1_IRGs_naming.csv | cut -d ',' -f6 > tempIRGprot.txt;
cat tempIRGprot.txt | sed '/^$/d' > temp2;
cat temp2 | sed '/?/d' > tempIRGprot.txt;
sed 's/\;/\n/g' tempIRGprot.txt > temp;
sed 's/ //g' temp > tempIRGprot.txt;
rm temp; rm temp2

# Extract all orthologues: 

head -n 1 $TOR_PATH/big_data/proteomes/Results_OrthoFinder/Phylogenetic_Hierarchical_Orthogroups/N0.tsv > $TOR_PATH/GIT/RodentIRGs/rawdata/orthoFinder_sub/N0_HOG_IRG.tsv;
rg -f $TOR_PATH/GIT/RodentIRGs/tempIRGprot.txt $TOR_PATH/big_data/proteomes/Results_OrthoFinder/Phylogenetic_Hierarchical_Orthogroups/N0.tsv >> $TOR_PATH/GIT/RodentIRGs/rawdata/orthoFinder_sub/N0_HOG_IRG.tsv
rm tempIRGprot
