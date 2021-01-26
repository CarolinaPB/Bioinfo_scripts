source /home/WUR/moiti001/miniconda3/etc/profile.d/conda.sh
conda activate align-genomes

for i in {0..30}
do
    

    echo $scaff
    mummerplot -r chr$i -p redundans_oldturkey_chr$i --png redundans_oldturkey.filter

done