#!/usr/bin/env python
import sys
import os
import re
from decimal import Decimal
import numpy as np

def extract(fqstat):
    q30 = []
    basenum = []
    readnum = []
    gccontent = []
    with open(fqstat,'r') as fin:
        for line in fin:
            line = line.strip('\n')
            searchObj1 = re.search(r'^#BaseQ:20--30.*>Q30:\s(.*)%$',line)
            searchObj2 = re.search(r'^#ReadNum:.*BaseNum:\s(\d+)',line)
            searchObj3 = re.search(r'^#ReadNum:\s(\d+)',line)
            searchObj4 = re.search(r'^#GC%:\s(\d+)',line)
            if searchObj1:
                q30.append(float(searchObj1.group(1)))
            if searchObj2:
                basenum.append(float(searchObj2.group(1)))
            if searchObj3:
                readnum.append(float(searchObj3.group(1)))
            if searchObj4:
                gccontent.append(float(searchObj4.group(1)))
    q30_mean = get_mean(q30)
    basenum_sum = np.sum(basenum)
    basenum_sum = Decimal(basenum_sum/1000000).quantize(Decimal("1."))
    readnum_sum = np.sum(readnum)
    gccontent_mean = get_mean(gccontent)
    return q30_mean,basenum_sum,readnum_sum,gccontent_mean
        
def get_mean(list):
    total = sum(list)
    length = len(list)
    mean  = Decimal(total/length).quantize(Decimal("0.00"))
    return mean

def map_reads(flagstat):
    total_num = 0 
#    mate_num = 0  
#    singleton_num = 0 
    mapped_num = 0 
    with open(flagstat,'r') as fin: 
        for line in fin: 
            line = line.strip('\n') 
            searchObj1 = re.search(r'(\d+)\s\+\s\d+\sin\stotal\s\(QC-passed\sreads\s\+\sQC-failed\sreads\)',line) 
#            searchObj2 = re.search(r'(\d+)\s\+\s\d+\swith\sitself\sand\smate\smapped',line) 
            searchObj3 = re.search(r'(\d+)\s\+\s\d+\smapped*',line) 
            if searchObj1: 
                total_num = float(searchObj1.group(1)) 
#            elif searchObj2: 
#                mate_num = float(searchObj2.group(1)) 
            elif searchObj3:  
                mapped_num = float(searchObj3.group(1)) 
            else: 
                continue 
#    total_map_num = int(mate_num + singleton_num) 
    ratio = Decimal(mapped_num*100/total_num).quantize(Decimal("0.00")) 
    return(ratio,int(mapped_num))

#tup = ('ITGB7_29','TBP_30','PUM1','POP4','BLOC1S2','RET_1',"RET_2")
#tup = ('CCDC6-RET_C1-R12','NCOA4-RET_N8-R12','ETV6-NTRK3_E4-N14','ITGB7_E5-E6','TBP_E6-E7','PUM1_E21-E22','POP4_E6-E7','BLOC1S2_E5-E6','KRT7_E4-E5','TG_E38-E39','SLC5A5_E14-E15')
#tup = ('ITGB7_E5-E6','TBP_E6-E7','PUM1_E21-E22','POP4_E6-E7','BLOC1S2_E5-E6','KRT7_E4-E5','TG_E38-E39','SLC5A5_E14-E15')
tup = ('ITGB7_E5-E6','TBP_E6-E7','PUM1_E21-E22','POP4_E6-E7','BLOC1S2_E5-E6')
def inter_ref_gene(reads_stat):
    cov_dic = {}
    with open(reads_stat,'r') as fin:
        for line in fin.readlines():
            line = line.strip('\n')
            gene = line.split('\t')[0]
            cov_depth = line.split('\t')[1]
            if gene in tup:
                cov_depth = str(Decimal(cov_depth).quantize(Decimal("1.")))
                cov_dic[gene] = cov_depth
            else:
                continue
#    imbalance = (float(cov_dic['RET_2']) - float(cov_dic['RET_1'])) / (float(cov_dic['ITGB7_29']) + float(cov_dic['TBP_30']) + float(cov_dic['PUM1']) + float(cov_dic['POP4']) + float(cov_dic['BLOC1S2']))
#    imbalance = Decimal(imbalance).quantize(Decimal("0.0000"))
    return (cov_dic)

def main(fqstat_raw,fqstat_clean,flagstat,reads_stat,samplename):
    q30_mean_raw,basenum_num_raw,readnum_sum_raw,gc_content_raw =  extract(fqstat_raw)
    q30_mean_clean,basenum_num_clean,readnum_sum_clean,gc_content_clean = extract(fqstat_clean)
    map_ratio,total_map_num = map_reads(flagstat)
    cov_dic = inter_ref_gene(reads_stat)
    ic_num = 0
    ic_read = 0
    cont_num = 0
    status = ''
    cov_list = []
    warning = []
    for gene in tup:
        if cov_dic.has_key(gene):
            if basenum_num_clean == 0:
                gene_num = int(cov_dic[gene])
            else:
                gene_num = Decimal(int(cov_dic[gene]) * 50 / basenum_num_clean).quantize(Decimal("1."))
            cov_list.append(str(gene_num))
            if int(gene_num) > 0:
                ic_num += 1
                ic_read += int(gene_num)
                if int(gene_num) >= 500:
                    cont_num = cont_num + 1
        else:
            cov_list.append('.')
    if cont_num < 2:
        warning.append("IC_num=%s"%cont_num)
    if q30_mean_raw < 80:
        warning.append("Raw_Q30=%s"%q30_mean_raw)
#    if readnum_sum_clean/1000000 < 0.3:
#        warning.append("Clean_Reads=%d"%readnum_sum_clean)
    if map_ratio < 90:
        warning.append("Map_Ratio=%s"%map_ratio)
    if len(warning) == 0:
        status = 'Pass'
    else:
        status = 'Fail'
#    print '\t'.join(samplename,basenum_num_raw,q30_mean_raw,basenum_num_clean,q30_mean_clean,total_map_num,map_ratio,cov_list)
#    print samplename,'\t',basenum_num_raw,'\t',readnum_sum_raw,'\t',basenum_num_clean,'\t',readnum_sum_clean,'\t',q30_mean_raw,'\t',q30_mean_clean,'\t',gc_content_clean,'\t',total_map_num,'\t',map_ratio,'\t','\t'.join(cov_list)

#    print '%s\t%d\t%.1f\t%d\t%.1f\t%.2f\t%.2f\t%.2f\t%d\t%.2f\t'%(samplename,basenum_num_raw,readnum_sum_raw,basenum_num_clean,readnum_sum_clean,q30_mean_raw,q30_mean_clean,gc_content_clean,total_map_num,map_ratio),'\t'.join(cov_list),"\t%s\t"%ic_read,';'.join(warning),"\t%s"%status
    print '%s\t%d\t%d\t%d\t%d\t%.2f\t%.2f\t%.2f\t%d\t%.2f\t'%(samplename,basenum_num_raw,readnum_sum_raw,basenum_num_clean,readnum_sum_clean,q30_mean_raw,q30_mean_clean,gc_content_clean,total_map_num,map_ratio),'\t'.join(cov_list),"\t%s\t"%ic_read,';'.join(warning),"\t%s"%status
if __name__ == '__main__':
    if len(sys.argv) - 1 < 4:
        print 'py fqstat_raw fqstat_clean flagstat reads_stat samplename '
        exit()
#    print 'SampleName\tRaw_Data\tQ30_Raw\tClean_Data\tQ30_Clean\tMap_Num\tMap_ratio\tCCDC6-RET\tNCOA4-RET\tETV6-NTRK3\tITGB7\tTBP\tPUM1\tPOP4\tBLOC1S2\tKRT7\tTG\tSLC5A5'
    print 'Sample_Name\tRaw_Data(Mb)\tRaw_Reads\tClean_Data(Mb)\tClean_Reads\tRaw_Q30\tClean_Q30\tGC_Content\tMapped_Num\tMapped_Ratio\tITGB7\tTBP\tPUM1\tPOP4\tBLOC1S2\tIC_total_reads\tWarning\tStatus'

    fqstat_raw = os.path.abspath(sys.argv[1])
    fqstat_clean = os.path.abspath(sys.argv[2])
    flagstat = os.path.abspath(sys.argv[3])
    reads_stat = os.path.abspath(sys.argv[4])
    samplename = sys.argv[5]
    main(fqstat_raw,fqstat_clean,flagstat,reads_stat,samplename)
