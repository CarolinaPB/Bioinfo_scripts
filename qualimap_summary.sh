#for single mapping
for i in "$@"
do
case $i in
    -p=*|--prefix=*)
    PREFIX="${i#*=}"

    ;;
    -f=*|--file=*)
    FILE="${i#*=}"
    ;;
esac
done

if [ "$PREFIX" == "" ]
    then
        PREFIX="file"
fi

if [ "$FILE" == "" ]
    then
        FILE="genome_results.txt"
fi


if [ ! -f $FILE ]; then
    echo "ERROR: $FILE not found!"
    exit
fi

RES=$(grep "bam" $FILE | awk -F '=' '{print $NF}')
mean_map_qual=$(grep "mean mapping quality" $FILE | awk -F '=' '{print $NF}')
coverage=$(grep "mean coverage" $FILE | awk -F '=' '{print $NF}')
stdcoverage=$(grep "std coverage" $FILE | awk -F '=' '{print $NF}')
map_rate=$(grep "number of mapped reads" $FILE | awk -F '=' '{print $NF}')
num_reads=$(grep "number of reads" $FILE | awk -F '=' '{print $NF}')
num_supp=$(grep "number of supplementary alignments" $FILE | awk -F '=' '{print $NF}')
num_secondary=$(grep "number of secondary alignments" $FILE | awk -F '=' '{print $NF}')
num_duplicated=$(grep "number of duplicated reads" $FILE | awk -F '=' '{print $NF}')
duplication_rate=$(grep "duplication rate" $FILE | awk -F '=' '{print $NF}')
mapping_quality=$(grep "mean mapping quality" $FILE | awk -F '=' '{print $NF}')
gc_percentage=$(grep "GC percentage" $FILE | awk -F '=' '{print $NF}')

 
printf 'File\t%s\n' $RES >> ${PREFIX}_mapquality_summary.tsv
printf 'Mean coverage\t%s\n' $coverage >> ${PREFIX}_mapquality_summary.tsv
printf 'Std coverage\t%s\n' $stdcoverage >> ${PREFIX}_mapquality_summary.tsv
printf 'Mapping quality\t%s\n' $mapping_quality >> ${PREFIX}_mapquality_summary.tsv
printf 'Number of reads\t%s\n' $num_reads >> ${PREFIX}_mapquality_summary.tsv
printf 'Number of mapped reads\t%s\n' "$map_rate" >> ${PREFIX}_mapquality_summary.tsv
printf 'Number of supplementary alignments\t%s\n' "$num_supp" >> ${PREFIX}_mapquality_summary.tsv
printf 'Number of secondary alignments\t%s\n' $num_secondary >> ${PREFIX}_mapquality_summary.tsv
printf 'Number of duplicated reads\t%s\n' $num_duplicated >> ${PREFIX}_mapquality_summary.tsv
printf 'Duplication rate\t%s\n' $duplication_rate >> ${PREFIX}_mapquality_summary.tsv
printf 'GC\t%s\n' $gc_percentage >> ${PREFIX}_mapquality_summary.tsv
