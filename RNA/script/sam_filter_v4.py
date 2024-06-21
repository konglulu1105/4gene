#!/usr/bin/env python
import sys
import os
import re


if len(sys.argv) < 3:
        print 'py  fusion.bed sample_sort.sam '
        exit(1)
bed_input = os.path.abspath(sys.argv[1])
sam_input = os.path.abspath(sys.argv[2])
#output_path = os.path.dirname(sam_input)

def fusion_bed(bed_input):
    dic_bed = {}
    arr = []
    with open(bed_input,'r') as fin:
        for line in fin:
            if line.startswith('#'):
                continue
            else:
                line = line.strip('\n')
                bed_key = str(line.split('\t')[0])
                arr.append(bed_key)
                bed_start = int(line.split('\t')[7])
                bed_end = int(line.split('\t')[8])
                insert_range = int(line.split('\t')[3])
                bed_break = int(line.split('\t')[6])
                dic_bed[bed_key] = (bed_start,bed_end,insert_range,bed_break)
    return (dic_bed,arr)

def sum_list(items):  
    sum_numbers = 0  
    for x in items:  
        sum_numbers += x  
    return sum_numbers

def bwa_sam(sam_input,bed_input):
    dic_bed = fusion_bed(bed_input)[0]
    dic_num = {}
#    readid_list = []
    sam_file = open(sam_input + '.filter','w')
    with open(sam_input,'r') as fin:
        for line in fin:
            line = line.strip('\n')
            if line.startswith('@'):
#                print line
                sam_file.write(line + '\n')
            else:
#                readid = line.split('\t')[0]
                rname = line.split('\t')[2]
                if rname not in dic_num.keys():
                    dic_num[rname] = 0
                pos = int(line.split('\t')[3])
                cigar = line.split('\t')[5]
                rnext = line.split('\t')[6]
                isize = abs(int(line.split('\t')[8]))
                searchObj_1 = re.findall(r'(\d{1,3})I',cigar)
                searchObj_2 = re.findall(r'(\d{1,3})M',cigar)
                searchObj_3 = re.findall(r'(\d{1,3})D',cigar)
                searchObj_4 = re.findall(r'(\d{1,3})S',cigar)
                total_insert = 0
                total_match = 0
                total_delete = 0
                total_soft = 0
                if searchObj_1:
                    searchObj_1 = [int(i) for i in searchObj_1]
                    total_insert = sum_list(searchObj_1)
                if searchObj_2:
                    searchObj_2 = [int(i) for i in searchObj_2]
                    total_match = sum_list(searchObj_2)
                if searchObj_3:
                    searchObj_3 = [int(i) for i in searchObj_3]
                    total_delete = sum_list(searchObj_3)
                if searchObj_4:
                    searchObj_4 = [int(i) for i in searchObj_4]
                    total_soft = sum_list(searchObj_4)
                if rname in dic_bed.keys():
                    (start,end,insert,break_point) = dic_bed[rname]
                    if total_match >= 80 and total_insert <= 5 and total_delete <=5 and total_soft <= 50 and rnext == '=' and \
                    dic_bed[rname][3] - pos >= 10 and pos + total_match + total_insert - dic_bed[rname][3] >= 10 and  \
                    abs(isize - dic_bed[rname][2]) <= 20 and total_match >= (dic_bed[rname][2] * 0.6):
                        dic_num[rname] =  dic_num[rname] + 1
                        sam_file.write(line + '\n')

    for key in dic_bed.keys():
        if key not in dic_num.keys():
            dic_num[key] = 0
    sam_file.close()
    return(dic_num) 

dic_num = bwa_sam(sam_input,bed_input)
arr_bed = fusion_bed(bed_input)[1]
for i in arr_bed:
    print i + '\t' + str(dic_num[i])
