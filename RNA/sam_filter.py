#!/usr/bin/env python3
import os
import sys
import argparse
from decimal import Decimal
import re

REAL_DIR = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(os.path.dirname(REAL_DIR), "common"))
import tools


def parse_arguments():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("--fusion-bed", dest="fusion_bed", required=True, help="fusion bed")
    parser.add_argument("--input-sam", dest="input_sam", required=True, help="input bam")
    parser.add_argument("--filter-sam", dest="filter_sam", required=True, help="filter sam")
    parser.add_argument("--out-stat", dest="out_stat",
                        required=True, help="output reads stat")
    return parser.parse_args()

def main():
    args = parse_arguments()
    bed_input = os.path.abspath(args.fusion_bed)
    input_sam = os.path.abspath(args.input_sam)
    filter_sam = os.path.abspath(args.filter_sam)
    out_stat = os.path.abspath(args.out_stat)
    out_reads_stat(bed_input, input_sam, filter_sam, out_stat)
    
    
def out_reads_stat(bed_input, input_sam, filter_sam, out_stat):
    fout = open(out_stat, 'w')
    arr_list, bed_dict = read_fusion_bed(bed_input)
    num_dict = bwa_sam(input_sam, filter_sam, bed_dict)
    for i in arr_list:
        fout.write(i + '\t' + str(num_dict[i]) + "\n")
    fout.close()
    

def read_fusion_bed(bed_input):
    """read fusion_ZC_new_v1.bed, which is located in RNA's database directory
    which has 9 columns, #Fusion_type       seqence_start  seqence_end  pcr_len  pcr_len_L  pcr_len_R  break_piont  primer_start  primer_end
    
    Args:
        bed_input (string): path to the fusion_ZC_new_v1.bed file

    Returns:
        bed_dict (List[string: tuple]): the dictionary of bed entries
        arr_list: list of bed keys
    """
    # # Convert TSV file to a list of dictionaries
    headers, entries = tools.tsv2list(bed_input)
    bed_dict = {}
    arr_list = []
    for entry in entries:
        bed_key = entry["#Fusion_type"]
        insert_range = int(entry["pcr_len"])
        bed_break = int(entry["break_piont"])
        bed_start = int(entry["primer_start"])
        bed_end = int(entry["primer_end"])
        arr_list.append(bed_key)
        # Create a tuple of values and assign it to the bed key in the dictionary
        bed_dict[bed_key] = (bed_start, bed_end,
                            insert_range, bed_break)
    return (arr_list, bed_dict)

def sum_list(items):
    sum_numbers = 0
    for x in items:
        sum_numbers += x
    return sum_numbers

def bwa_sam(input_sam, filter_sam, bed_dict):
    num_dict = {}
    sam_file = open(filter_sam, 'w')
    with open(input_sam, 'r') as fin:
        for line in fin:
            line = line.strip('\n')
            if line.startswith('@'):
                #                print line
                sam_file.write(line + '\n')
            else:
                #                readid = line.split('\t')[0]
                rname = line.split('\t')[2]
                if rname not in num_dict.keys():
                    num_dict[rname] = 0
                pos = int(line.split('\t')[3])
                cigar = line.split('\t')[5]
                rnext = line.split('\t')[6]
                isize = abs(int(line.split('\t')[8]))
                search_obj1 = re.findall(r'(\d{1,3})I', cigar)
                search_obj2 = re.findall(r'(\d{1,3})M', cigar)
                search_obj3 = re.findall(r'(\d{1,3})D', cigar)
                search_obj4 = re.findall(r'(\d{1,3})S', cigar)
                total_insert = 0
                total_match = 0
                total_delete = 0
                total_soft = 0
                if search_obj1:
                    search_obj1 = [int(i) for i in search_obj1]
                    total_insert = sum_list(search_obj1)
                if search_obj2:
                    search_obj2 = [int(i) for i in search_obj2]
                    total_match = sum_list(search_obj2)
                if search_obj3:
                    search_obj3 = [int(i) for i in search_obj3]
                    total_delete = sum_list(search_obj3)
                if search_obj4:
                    search_obj4 = [int(i) for i in search_obj4]
                    total_soft = sum_list(search_obj4)
                if rname in bed_dict:
                    (start, end, insert, break_point) = bed_dict[rname]
                    if (total_match >= 80 and total_insert <= 5 and total_delete <= 5 
                        and total_soft <= 50 and rnext == '=' and bed_dict[rname][3] - pos >= 10 
                        and pos + total_match + total_insert - bed_dict[rname][3] >= 10 and abs(isize - bed_dict[rname][2]) <= 20 
                        and total_match >= (bed_dict[rname][2] * 0.6)):
                        num_dict[rname] = num_dict[rname] + 1
                        sam_file.write(line + '\n')
    for key in bed_dict.keys():
        if key not in num_dict:
            num_dict[key] = 0
    sam_file.close()
    return num_dict

if __name__ == '__main__':
    main()