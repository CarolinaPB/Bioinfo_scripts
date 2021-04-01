#plot SV larger than minsize from a SV vcf summary table

# vcf-query -f '%CHROM\t%POS\t%INFO/END\t%ID\t%INFO/SVLEN\n' F1/variants.vcf | sed "s/svim.//g" | awk '{sub(/\..*$/,"", $4)}1' OFS=\\t  | sed "s/\-//g" | sort -n -k5,5 > variant_table.txt
# awk '$1 ~ /^[[:digit:]]+$/'
re='^[0-9]+$'


minsize=900000
cat variant_table.txt |  while read line
do
    chr=$(echo $line | awk '{print $1}')
    start=$(echo $line | awk '{print $2}')
    end=$(echo $line | awk '{print $3}')
    type=$(echo $line | awk '{print $4}')
    size=$(echo $line | awk '{print $5}')
    if  [[ $size =~ $re ]] ; then
        if [ $size -gt  $minsize ]; then
            echo $chr $start $end $type
            samplot plot \
            -n F1 parent1 parent2 \
            -b /lustre/nobackup/WUR/ABGC/moiti001/results/Mgal_WUR_HG_1.0/structural_var/F1_15022021/map_ngmlr_F1.sorted.bam \
            /lustre/nobackup/WUR/ABGC/moiti001/results/Mgal_WUR_HG_1.0/structural_var/parent1/map_ngmlr_parent1.sorted.bam \
            /lustre/nobackup/WUR/ABGC/moiti001/results/Mgal_WUR_HG_1.0/structural_var/parent2/map_ngmlr_parent2.sorted.bam \
            -c $chr \
            -t $type \
            --start $start \
            --end $end \
            --same_yaxis_scales 
        fi
    fi
done