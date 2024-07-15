#!/usr/bin/env python3
import os
import sys
import argparse
import re
REAL_DIR = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(os.path.dirname(REAL_DIR), "common"))
import tools


def parse_arguments():
    parser = argparse.ArgumentParser(description="combine snv vcf with indel vcf file")
    parser.add_argument("--snv-vcf", dest="snv_vcf", required=True,
                        help="vcf format file that contains snvs")
    parser.add_argument("--indel-vcf", dest="indel_vcf", required=True,
                        help="vcf format file that contains indels")
    parser.add_argument("--out-vcf", dest="out_vcf", required=True, help="")
    return parser.parse_args()

def main():
    args = parse_arguments()
    snv_vcf = args.snv_vcf
    indel_vcf = args.indel_vcf
    out_vcf = args.out_vcf
    combine_vcf(snv_vcf, indel_vcf, out_vcf)
    
def combine_vcf(snv_vcf, indel_vcf, out_vcf):
    meta_info, snps = load_vcf(snv_vcf)
    meta_info, indels = load_vcf(indel_vcf)
    fout = open(out_vcf, 'w')
    fout.write("\n".join(meta_info)+"\n")
    fout.write("\n".join(snps + indels) + "\n")
    fout.close()

def load_vcf(vcf_file):
    meta_info = list()
    content = list()
    var_dict = dict()
    with open(vcf_file, 'r') as finp:
        for line in finp.readlines():
            line = line.strip("\n")
            if re.match(r'^#', line):
                meta_info.append(line)
                continue
            if re.match(r'^chr', line):
                cols = line.split("\t")
                key = (cols[0], cols[1], cols[3], cols[4])
                var_dict[key] = line
                content.append(line)
    return meta_info, content
        
    
    
if __name__ == '__main__':
    main()