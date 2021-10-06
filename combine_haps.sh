out="combined.fasta"
cat chr5_393462_677667.fasta > $out

for dir in polish_output/*; do
   pid=`basename $dir`

   echo ">${pid}.hapA.fasta" >> $out
   cat ${dir}/*hapA.polished.fasta | tail -n +2 >> $out

   echo ">${pid}.hapB.fasta" >> $out
   cat ${dir}/*hapB.polished.fasta | tail -n +2 >> $out

   echo $pid
done
