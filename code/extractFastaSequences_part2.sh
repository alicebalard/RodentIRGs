#!/bin/bash

echo "Hello! which table do we chose? First column must be user names (e.g. Mmus_Irgb2), second NCBI names"
read -e -p "My table:" mytable

echo "How will we call our output fasta file?"
read namefastaoutput

echo "We look for sequences in NCBI, using esearch and efetch from Entrez NCBI package:"
cut -d"," -f2 $mytable > temp1

while read line; 
	do esearch -db protein -query $line < /dev/null | efetch -format fasta >> $namefastaoutput; done < temp1
echo "DONE"

# Now, we make the output a one-liner
awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' $namefastaoutput > temp2
mv temp2 $namefastaoutput

# Prepare add-on for header
cut -d"," -f1 $mytable > addonHeader

#echo "Get one every two line starting from 1st:headers"
awk 'NR%2==1' $namefastaoutput > headers

# remove >
sed 's/>//' headers > headers2

# merge to make full headers
paste -d " " addonHeader headers2 > fullheaders

# add > at every line start
sed 's/^/>/' fullheaders > fullHeaders

# one every two line starting from 0th: sequences
awk 'NR%2==0' $namefastaoutput > sequences

# merge line by line using headers from the text file
paste -d'\n' fullHeaders sequences > $namefastaoutput

#rm trash
rm temp1; rm fullheaders; rm fullHeaders; rm sequences; rm headers; rm headers2 ; rm addonHeader
