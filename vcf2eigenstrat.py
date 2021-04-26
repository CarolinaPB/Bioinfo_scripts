# Convert a vcf file to eigenstrat format
# removes multi-alleleic and indel sites. 
# usage: python vcf2eigenstrat.py -v vcf_file.vcf(.gz) -o out_root
# will generate out_root.[snp,ind,geno].
# removed multiallelic sites and indels
# Deals with haploid cases including mixed haploid/diplod like X as well. 
# -i option is a .ind file to get population names and sex. 

# adapted from https://github.com/mathii/gdc/blob/master/vcf2eigenstrat.py

from __future__ import division
import sys, vcf
import argparse, gzip

################################################################################

def parse_options():
    """
    Options are described by the help() function
    """
    options ={ "vcf_file":None, "out":"out", "ref":None, "indAsPop":False, "indmap":None  }
    parser = argparse.ArgumentParser( description='Convert vcf to eigenstrat formats')
    parser.add_argument("-v", "--vcf_file", help="VCF_file (compressed)", nargs=1, default=None)
    parser.add_argument("-o", "--out", help="outfile prefix", nargs=1, default="out")
    parser.add_argument("-ref", "--ref", help="ref", nargs=1, default=None)
    parser.add_argument("-indAsPop", "--indAsPop", help="indAsPop", nargs=1, default=None)
    parser.add_argument("-indmap", "--indmap", help="indmap", nargs=1, default=None)
    
    args = parser.parse_args()
    
    options["vcf_file"] = args.vcf_file[0]
    options["ref"] = args.ref
    options["indmap"] = args.indmap
    options["indAsPop"] = args.indAsPop
    options["out"] = args.out[0]
    
    print("found options:")
    print(options)
    return(options)


################################################################################

def main(options):
    """
    Convert vcf to eigenstrat format (ind, snp and geno files)
    """
    # vcf=gdc.open2(options["vcf"])
    # vcf= vcf.Reader(filename = options["vcf_file"], strict_whitespace=True)
    with gzip.open(options["vcf_file"], "rb") as vcf:
        # snp = open(options["out"] + ".snp")
        # ind = open(options["out"] + ".ind")
        # geno = open(options["out"] + ".geno")
        snp, ind, geno = [open(options["out"]+x, "w") for x in [".snp", ".ind", ".geno"]]
        removed={"multiallelic":0, "indel":0}
        count=0
        
        if options["indmap"]:
            pop_map={}
            sex_map={}
            ind_map_file=open(options["indmap"], "r")
            for line in ind_map_file:
                bits=line[:-1].split()
                pop_map[bits[0]]=bits[2]
                sex_map[bits[0]]=bits[1]
            ind_map_file.close()
        
        for line in vcf:
            if line.startswith(b"##"):				  # Comment line
                next
            elif line[:6]==b"#CHROM":			  # Header line
                # print(line)
                inds=line.split()[9:]
                if options["ref"]:
                    # print("here-1")
                    ind.write(options["ref"]+"\tU\tREF\n")
                
                if options["indmap"]:
                    for indi in inds:
                        # print("here0")
                        ind.write(indi.decode("utf-8") +"\t"+sex_map.get(indi, "U")+"\t"+pop_map.get(indi, "POP")+"\n")
                elif options["indAsPop"]:
                    for indi in inds:
                        # print("here1")
                        ind.write(indi.decode("utf-8")+"\tU\t"+indi.decode("utf-8") +"\n")
                else:
                    for indi in inds:
                        # print("here2")
                        # print(indi.decode("utf-8") +"\tU\tPOP\n")
                        ind.write(indi.decode("utf-8") +"\tU\tPOP\n")
                    
            else:							  # data
                bits=line.split()

                if b"," in bits[4]:
                    removed["indel"]+=1
                    continue
                if len(bits[3])!=1 or len(bits[4])!=1:
                    removed["multiallelic"]+=1
                    continue
                else:
                    if bits[2]==b".":
                        bits[2]=bits[0]+b":"+bits[1]
                    snp.write("    ".join([bits[2].decode("utf-8"), bits[0].decode("utf-8") , "0.0", bits[1].decode("utf-8") , bits[3].decode("utf-8") , bits[4].decode("utf-8") ])+"\n")
                    geno_string=""
                    if options["ref"]:
                        geno_string="2"
                    for gt in bits[9:]:
                        geno_string+=decode_gt_string(gt)
                    geno.write(geno_string+"\n")
                    count+=1

        [f.close for f in [ind, snp, geno]]

        print("Done. Wrote "+str(count) + " sites")
        print("Excluded " + str(sum(removed.values())) + " sites")
        for key in removed:
            print("Excluded " + str(removed[key]) + " " + key)
        return

################################################################################

def decode_gt_string(gt_string):
    """
    Tries to work out the genotype from a vcf genotype entry. 9 for missing [or not in {0,1,2}]
    """
    gt=gt_string.split(b":")[0]
    gt = gt.decode("utf-8") 
    # print(gt)
    # print(len(gt))
    
    if len(gt)==1:
        if gt=="0":                       # haploid
            return "2"
        elif gt=="1":
            return "0"
        else:
            return "9"
    elif len(gt)==3:
        if gt[0]=="0" and gt[2]=="0":
            return "2"
        if gt[0]=="0" and gt[2]=="1":
            return "1"
        if gt[0]=="1" and gt[2]=="0":
            return "1"
        if gt[0]=="1" and gt[2]=="1":
            return "0"
        else:
            return "9"

    raise Exception("Unknown genotype: "+gt)
        
################################################################################

if __name__=="__main__":
	options=parse_options()
	main(options)
	