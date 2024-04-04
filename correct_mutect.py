#!/usr/bin/python3

# from: https://github.com/broadinstitute/gatk/issues/6857#issuecomment-1246662957

# makes Mutect2 output vcf compliant with VCF specification in case of AS_FilterStatus delimiters at multiallelic loci
# ("wrong number of fields in AS_FilterStatus?")

# https://github.com/broadinstitute/gatk/issues/6857

# usage: cat sample.vcf | correct_mutect.py | bgzip > sample.right.vcf.gz && tabix -p vcf sample.right.vcf.gz

import sys

if __name__ == "__main__":
    for l in sys.stdin:
        if l[0] != "#":
            cols = l.split("\t")
            # if "," in cols[4]: # not sure if the correction is needed at biallelic sites
            if True:
                info = cols[7].split(";")
                for i in range(len(info)):
                    a = info[i]
                    if a[:16] == "AS_FilterStatus=":
                        b = list(map(lambda c: c.split(","), a[16:].split("|")))
                        info[i] = "AS_FilterStatus=" + ",".join(map(lambda c: "|".join(c), b))
                cols = cols[:7] + [";".join(info)] + cols[8:]
                l = "\t".join(cols)
        print(l, end = "")
