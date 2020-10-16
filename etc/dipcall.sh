
# dipcall downloaded from here
# https://github.com/lh3/dipcall

dipcall="/home/scott/bin/dipcall.kit/run-dipcall"

ref=chr5_393462_677667.fasta

mkdir -p dipcall

function format {

   grep '#' $1 > __head__.txt
   grep -v '#' $1 > __body__.txt
   
   awk -v FS='\t' -v OFS='\t' '$2+=393462' __body__.txt | \
      sed 's/chr5_393462_677667/chr5/g'> __body2__.txt
   
   cat __head__.txt __body2__.txt > $2
   rm __head__.txt
   rm __body*__.txt
}

for dir in polish_output/*; do
   iid=`basename $dir`
   echo $iid
   mkdir -p dipcall/${iid}
   
   prefix=dipcall/${iid}/${iid}
   
   A=chr5_393462_677667hapA.polished.fasta
   B=chr5_393462_677667hapB.polished.fasta
   
   $dipcall $prefix $ref ${dir}/${A} ${dir}/${B} > ${prefix}.mak
   sed -i 's/xasm5/xasm20/g' ${prefix}.mak
   
   make -j2 -f ${prefix}.mak
   
   gunzip ${prefix}.pair.vcf.gz
   format ${prefix}.pair.vcf  ${prefix}.pair.vcf2
   mv ${prefix}.pair.vcf2 ${prefix}.pair.vcf
   gzip ${prefix}.pair.vcf
   rm ${prefix}.dip.vcf.gz
   make -j2 -f ${prefix}.mak

   gunzip ${prefix}.dip.vcf.gz
   sed -i "s/syndip/$iid/g" ${prefix}.dip.vcf
   bgzip ${prefix}.dip.vcf
   tabix ${prefix}.dip.vcf.gz

done
