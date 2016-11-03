#This script parses a taxon's BLAST table for its coresponding rDNA hits, and pulls them out of the assembly and orients the sequence correctly

cat TaxonList | while read taxon; do

        awk '{print $1"\t"$2}' ${taxon}_all_xanth_prot > rDNAContigNames
                cat rDNAContigNames | while read contig; do
                        bestContigHit=($(grep "$contig" ${taxon}_all_xanth_prot | head -n 1))
#                       echo ${bestContigHit[@]}
                        echo "$contig" > tmp3
                        actualContig=`awk '{print $1}' tmp3`
                        compressedContig=`sed 's/\t//g' tmp3`
                        hitBeg=${bestContigHit[6]}
                        hitEnd=${bestContigHit[7]}
                        diff=$((hitEnd - hitBeg))
#                       echo $hitBeg $hitEnd $diff
                        #Reverse complement sequence if necessary
                        if [ $diff -lt 0 ]; then
                                #get minumum and maximum positions of contig
                                minpos=`grep "$contig" ${taxon}* | awk '{print $8}' | sort -nk 1 | head -n 1`
                                maxpos=`grep "$contig" ${taxon}* | awk '{print $7}' | sort -nk 1 | tail -n 1`
                                grep -A 1 "$actualContig" /home/kyle/LichenWGS/AllAssembledGenomes/${taxon}_assembly.fasta > tmp
                                cat tmp | tail -n 1 | sed -e 's/\(.\)/\1 /g' | cut -f $minpos-$maxpos -d " " | sed 's/ //g' > forwardSeq
                                echo "Reverse-complemented ${taxon} hit ${contig}"
                                perl rc.pl `cat forwardSeq` > reversedSeq
                                echo "" >> reversedSeq
                                echo ">${taxon}_${compressedContig}" > tmp2
                                cat tmp2 reversedSeq > ${taxon}_${compressedContig}.fa
                                rm tmp; rm tmp2; rm reversedSeq; rm forwardSeq
                        else
                                #get minumum and maximum positions of contig
                                minpos=`grep "$contig" ${taxon}* | awk '{print $7}' | sort -nk 1 | head -n 1`
                                maxpos=`grep "$contig" ${taxon}* | awk '{print $8}' | sort -nk 1 | tail -n 1`
                                grep -A 1 "$actualContig" /home/kyle/LichenWGS/AllAssembledGenomes/${taxon}_assembly.fasta | tail -n 1 > tmp
                                cat tmp | tail -n 1 | sed -e 's/\(.\)/\1 /g' | cut -f $minpos-$maxpos -d " " | sed 's/ //g' > forwardSeq
                                echo ">${taxon}_${compressedContig}" > tmp2
                                cat tmp2 forwardSeq > ${taxon}_${compressedContig}.fa
                                rm tmp; rm tmp2; rm forwardSeq
                        fi
                        rm tmp3
                done
        rm rDNAContigNames
done
