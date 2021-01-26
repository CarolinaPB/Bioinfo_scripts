
cat /lustre/nobackup/WUR/ABGC/moiti001/results/scaffold_long_reads_lrscaf/scaff-salsa/lrscaf_top40_lengths.txt |  while read line
do
    scaff=$(echo $line | awk '{print $1}')

    echo $scaff
    mummerplot -r $scaff -p LRSCAF_chicken_$scaff --png LRSCAF_chicken.filter

done