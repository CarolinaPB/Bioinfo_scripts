module load bcftools

# usage: sh get_sniffles_summary.sh -vcf=VCF -o=OUTFILE


for i in "$@"
do
case $i in
    -vcf=*|--vcf=*)
    VCF="${i#*=}"
    ;;
esac
done


printf 'Chr\tSV\tStart\tEnd\tEnd_chr\tLen\tRead_support\tGT\n'> sniffles_summary.txt
bcftools query -f '%CHROM\t%INFO/SVTYPE\t%POS\t%INFO/END\t%INFO/CHR2\t%INFO/SVLEN\t%INFO/RE[\t%GT]\n' $VCF  >> sniffles_summary.txt
