import sys
if len(sys.argv) < 4:
        print 'py  fusion_reads.stats fusion_anno samplename '
        exit(1)

dic1 = {}

for i in open(sys.argv[2]).readlines():
    i = i.strip('\n')
    fusion_id = i.split('\t')[0]
    fusion_gene = i.split('\t')[1]
    break_point1 = i.split('\t')[2]
    break_point2 = i.split('\t')[3]
    tran1 = i.split('\t')[4]
    tran2 = i.split('\t')[5]
    dic1[fusion_id] = [fusion_gene,tran1,tran2,break_point1,break_point2]
print 'Sample_Name\tGene\tTranscript_1\tTranscript_2\tExon_Rank_1\tExon_Rank_2\tSupport_Reads\tStatus'
#for i in open(sys.argv[1]).readlines()[1:]:
for i in open(sys.argv[1]).readlines():
    i = i.strip('\n')
#    fusion_type = i.split('\t')[0].split(':')[0]
    fusion_id = i.split('\t')[0]
    reads_num = int(i.split('\t')[1])
#    cov_depth = int(cov_depth)
    if reads_num >= 8 and dic1.has_key(fusion_id):
        print '%s\t%s\t%s\t%s\t%s\t%s\t%d\t%s' %(sys.argv[3],dic1[fusion_id][0],dic1[fusion_id][1],dic1[fusion_id][2],dic1[fusion_id][3],dic1[fusion_id][4],reads_num,'NA')
