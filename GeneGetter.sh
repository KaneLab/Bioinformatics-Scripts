cat GenesInAll10 | while read gene; do

>${gene}.tmp

        cat directories | while read folder; do
                grep -h -A 1 $gene ${folder}/Single* > temp
                sed "s/$gene/$folder/g" temp >> ${gene}.tmp
                sed 's/\-\-//' ${gene}.tmp | sed '/^\s*$/d' > temp2; mv temp2 ${gene}
        done

        perl translatorX.pl -t T -i ${gene} -o ${gene}Aligned
done

rm *.tmp
rm *.muscle.log
rm *Aligned.html
rm *Aligned.nt12_ali.fasta
rm *Aligned.nt1_ali.fasta
rm *Aligned.nt2_ali.fasta
rm *Aligned.nt3_ali.fasta
#rm *Aligned.nt_ali.fasta
rm *Aligned.aaseqs
rm *Aligned.aa_based_codon_coloured.html
rm *Aligned.aaseqs.fasta

for i in `ls *nt_ali.fasta`; do

        bash fastaUnwrap.sh $i

done

rm *nt_ali.fasta
rm *aa_ali.fasta

for j in `ls *unwrapped.fa`; do

        lines=`wc -l $j | awk '{print $1}'`
        if [ $lines -ne 19 ]; then
                :
        else
                grep ">" $j > headerLine
                grep -v ">" $j > sequenceLine
                paste headerLine sequenceLine | sort -k1,1 | sed 's/    /\n/g' > ${j}.sorted
        fi
done

rm *unwrapped.fa
>seq
for k in `ls *.sorted`; do

        paste $k seq > temp3; mv temp3 seq

done

##################################################
transpose(){
awk '
{

    for (i=1; i<=NF; i++)  {

       a[NR,i] = $i

    }

}

NF>p { p = NF }

END {

    for(j=1; j<=p; j++) {

       str=a[1,j]

       for(i=2; i<=NR; i++){

           str=str" "a[i,j];

       }

       print str

    }

}' $1
}
##################################################


#formatting to get loci for which each taxon has data
sort headerLine > headerLineSorted
grep -v ">" seq | sed 's/       //g' > seq2; mv seq2 seq
paste headerLineSorted seq > seq2; mv seq2 AllProteinsConcatenated.fasta
awk '{print $2}' AllProteinsConcatenated.fasta | sed 's/\(.\{1\}\)/\1      /g' > temp
transpose temp > temp2
grep -v "-" temp2 > temp
transpose temp > temp2
sed 's/ //g' temp2 > temp
paste headerLineSorted temp > temp2
sed 's/\t/\n/g' temp2 > AllProteinsConcatenated_ready.fasta
