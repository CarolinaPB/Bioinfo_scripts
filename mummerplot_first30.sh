#!/bin/bash

# command line arguments 
# usage: sh mummerplot_first30.sh -p=PREFIX -f=FILTERFILE


for i in "$@"
do
case $i in
    -p=*|--prefix=*)
    PREFIX="${i#*=}"

    ;;
    -f=*|--filter=*)
    FILTER="${i#*=}"
    ;;
esac
done

echo $PREFIX
echo $FILTER
source /home/WUR/moiti001/miniconda3/etc/profile.d/conda.sh
conda activate align-genomes



# for i in {1..30}
for i in {1..35}
do
    

    echo $scaff
    # mummerplot -q Chr$i -p ${PREFIX}_chr$i --png $FILTER
    # mummerplot -q $i -p ${PREFIX}_chr$i --png $FILTER
    mummerplot -r $i -p ${PREFIX}_chr$i --png $FILTER

done

# mummerplot -q Chr41 -p ${PREFIX_chr}_Z --png $FILTER
# mummerplot -q Z -p ${PREFIX}_chrZ --png $FILTER
mummerplot -r Z -p ${PREFIX}_chrZ --png $FILTER