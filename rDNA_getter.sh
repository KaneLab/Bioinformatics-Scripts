#This script parses a taxon's BLAST table for its coresponding rDNA hits, and pulls them out of the assembly and orients the sequence correctly

cat TaxonList | while read taxon; do

        awk '{print $1}' ${taxon}* | sort | uniq -c | sed 's/_/\t/g' | awk '$1 > 3 && $5 > 2000' | sed 's/\t/_/g' | awk '{print $2}' > rDNAContigNames
                cat rDNAContigNames | while read contig; do
                        bestContigHit=($(grep $contig ${taxon}_photo_hits | head -n 1))
#                       echo ${bestContigHit[@]}
                        hitBeg=${bestContigHit[8]}
                        hitEnd=${bestContigHit[9]}
                        diff=$((hitEnd - hitBeg))
#                       echo $hitBeg $hitEnd $diff
                        #Reverse complement sequence if necessary
                        if [ $diff -lt 0 ]; then
                                #get minumum and maximum positions of contig
                                minpos=`grep $contig ${taxon}* | awk '{print $7}' | sort -nk 1 | head -n 1`
                                maxpos=`grep $contig ${taxon}* | awk '{print $8}' | sort -nk 1 | tail -n 1`
                                grep -A 1 $contig ../${taxon} > tmp
                                cat tmp | tail -n 1 | sed -e 's/\(.\)/\1 /g' | cut -f $minpos-$maxpos -d " " | sed 's/ //g' > forwardSeq
                                echo "Reverse-complemented ${taxon} hit ${contig}"
                                perl rc.pl `cat forwardSeq` > reversedSeq
                                echo "" >> reversedSeq
                                echo ">${taxon}_${contig}" > tmp2
                                cat tmp2 reversedSeq > ${taxon}_${contig}.fa
                                rm tmp; rm tmp2; rm reversedSeq; rm forwardSeq
                        else
                                #get minumum and maximum positions of contig
                                minpos=`grep $contig ${taxon}* | awk '{print $7}' | sort -nk 1 | head -n 1`
                                maxpos=`grep $contig ${taxon}* | awk '{print $8}' | sort -nk 1 | tail -n 1`
                                grep -A 1 $contig ../${taxon} | tail -n 1 > tmp
                                cat tmp | tail -n 1 | sed -e 's/\(.\)/\1 /g' | cut -f $minpos-$maxpos -d " " | sed 's/ //g' > forwardSeq
                                echo ">${taxon}_${contig}" > tmp2
                                cat tmp2 forwardSeq > ${taxon}_${contig}.fa
                                rm tmp; rm tmp2; rm forwardSeq
                        fi
                done
        rm rDNAContigNames
done
