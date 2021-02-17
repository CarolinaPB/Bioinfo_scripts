printf 'Sample\tMean coverage\tMapping rate\tMean mapping quality\n'> sample_quality_summary.tsv

for dir in qualimap/*; do
    if [ -d "$dir" ]; then
        for sample in $dir; do
            # echo $sample
            if [ -d "$sample" ]; then

                sample_name=$(echo $sample | awk -F '/' '{print $NF}')

                mean_map_qual=$(grep "mean mapping quality" $sample/genome_results.txt | awk -F '=' '{print $NF}')
                coverage=$(grep "mean coverage" $sample/genome_results.txt | awk -F '=' '{print $NF}')
                map_rate=$(grep "number of mapped reads" $sample/genome_results.txt | awk -F '(' '{print $NF}' | grep -o '[0-9.%]\+')

                printf '%s\t%s\t%s\t%s\n' $sample_name $coverage $map_rate $mean_map_qual >> sample_quality_summary.tsv
 
            fi
        done
    fi
done