module load vcftools

# command line arguments 
# usage: sh get_struct_var_vcf_subset.sh -vcf=VCF -o=OUTFILE


for i in "$@"
do
case $i in
    -vcf=*|--vcf=*)
    VCF="${i#*=}"
    ;;
    -o=*|--outfile=*)
    OUTFILE="${i#*=}"
    ;;
    -c=*|--svcaller=*)
    SVCALLER="${i#*=}"
    ;;
esac
done
#"s/svim_asm.//g"

if [ $SVCALLER == "svim" ]; then
    toreplace="svim."
elif [ $SVCALLER == "svim-asm" ]; then
    toreplace="svim_asm."
fi

vcf-query -f '%CHROM\t%ID\t%INFO/SVLEN\n' $VCF | sed "s/$toreplace//g" | awk '{sub(/\..*$/,"", $2)}1' OFS=\\t > $OUTFILE
