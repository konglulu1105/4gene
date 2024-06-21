#!/usr/bin/env python
import sys

if len(sys.argv) < 2:
    print "py sample_vep.txt"
    exit()
tup = ('Chr','Position_Start','Position_End','Reference','Mutation','GT','DP','AD',
'Mutation_Frequency','Gene','HGVSc','HGVSp','Transcript','Canonical','Exon|Intron','Mutation_Type',
'dbSNP','Cosmic','1000G_AF','1000G_EAS_AF','gnomAD_AF','gnomAD_EAS_AF','SIFT_score','SIFT','PolyPhen_score',
'PolyPhen','CLIN_SIG','simpleRepeat','Dedup_GT','Dedup_DP','Dedup_AD','Dedup_AF')
dic = {'A':'T','T':'A','G':'C','C':'G'}
tup_1 = ('A','T','G','C')
with open(sys.argv[1],'r') as fin:
        for line in fin:
            line = line.strip('\n')
            if line.startswith('Chr'):
                print line
            else:
                chrom = line.split('\t')[tup.index('Chr')]
                start = int(line.split('\t')[tup.index('Position_Start')])
                end = int(line.split('\t')[tup.index('Position_End')])
                reference = line.split('\t')[tup.index('Reference')]
                mutation = line.split('\t')[tup.index('Mutation')]
                gene = line.split('\t')[tup.index('Gene')]
                transcript = line.split('\t')[tup.index('Transcript')]
                hgvsc = line.split('\t')[tup.index('HGVSc')]
                hgvsp = line.split('\t')[tup.index('HGVSp')]
                if reference in tup_1 and mutation in tup_1:
                    if gene == 'TERT' and  start == end and start > 1295104 and end < 1295286 and hgvsc == '.' and hgvsp == '.':
                        site = 1295104 - start
                        hgvsc = '{0}:c.{1}{2:s}>{3:s}'.format(transcript,site,dic[reference],dic[mutation])
                        print '\t'.join(line.split('\t')[0:tup.index('HGVSc')]) + '\t' + hgvsc + '\t' + '\t'.join(line.split('\t')[tup.index('HGVSp'):])
                    else:
                        print line
                else:
                    print line
