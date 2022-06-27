#! /bin/bash
# from https://gist.github.com/slavailn/cafd7c50276a37b3aad61756ba994353

# 22498713 reads; of these:
#   22498713 (100.00%) were unpaired; of these:
#     1754404 (7.80%) aligned 0 times
#     17667104 (78.52%) aligned exactly 1 time
#     3077205 (13.68%) aligned >1 times
# 92.20% overall alignment rate

echo -e "sampleID\ttotal_reads\tunmapped\taligned_one_time\taligned_multiple\talignment_rate"
for file in ./*.summary
do
    filename=`basename $file`
    samplename=${filename%.summary}
    
    total_reads=""
    unmapped=""
    aligned_one_time=""
    aligned_multiple=""
    alignment_rate=""
   
    while IFS= read -r line
    do
        if [[ $line == *reads* ]];
        then
            total_reads=`echo $line | awk '{print $1}'`;
        fi
        
        if [[ $line == *aligned[[:space:]]*0[[:space:]]*times* ]];
        then 
            unmapped_reads=`echo $line | awk '{print $1}'`;
        fi

        if [[ $line == *aligned[[:space:]]*exactly[[:space:]]*1[[:space:]]*time* ]]
        then
           aligned_one_time=`echo $line | awk '{print $1}'`;
        fi

        if [[ $line == *aligned[[:space:]]*\>1[[:space:]]*times*  ]]
        then
           aligned_multiple=`echo $line | awk '{print $1}'`;
        fi

        if [[ $line == *overall[[:space:]]*alignment[[:space:]]*rate*  ]]
        then
           alignment_rate=`echo $line | awk '{print $1}'`;
           alignment_rate=`echo $alignment_rate | sed s/%//`;
        fi
      
    done < "$file"
    echo -e "$samplename\t$total_reads\t$unmapped_reads\t$aligned_one_time\t$aligned_multiple\t$alignment_rate"
done
