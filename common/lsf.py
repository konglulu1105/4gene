#!/usr/bin/env python3
import os
import sys

def lsf_sub(script_path, sample_name, cpu_num):
    '''
    creat directory  with script name

    Args:
     script_path: path of script
     sample_name: name of sample
     cpu_num: num of cpu
    '''
    abs_path = os.path.abspath(script_path)
    script_file_name = os.path.basename(abs_path)
    script_file_dir = os.path.dirname(abs_path)
    bsublog_dir = os.path.join(script_file_dir,'lsf_' + script_file_name)
    my_command = 'mkdir -p {0:s}'.format(bsublog_dir)
    print(my_command)
    os.system(my_command)

    print(bsublog_dir)
    my_command = "export LSF_ENVDIR=/data01/lsf/conf; export LSF_SERVERDIR=/data01/lsf/10.1/linux3.10-glibc2.17-x86_64/etc; export LSF_LIBDIR=/data01/lsf/10.1/linux3.10-glibc2.17-x86_64/lib ;export LSF_BINDIR=/data01/lsf/10.1/linux3.10-glibc2.17-x86_64/bin; /data01/lsf/10.1/linux3.10-glibc2.17-x86_64/bin/bsub -n {0:s} -R span[ptile=24] -J {1:s} -q normal -m 'node01+2 node02+1 node03' -L /bin/bash -o {2:s}/{3:s}.out -e {2:s}/{3:s}.err sh {4:s}".format(cpu_num,sample_name,bsublog_dir,script_file_name,abs_path)
    print(my_command)
    os.system(my_command)
    
if __name__ == "__main__":
    if len(sys.argv) - 1 < 3:
        print('file, samplename, cpu')
        exit()
    abs_path = os.path.abspath(sys.argv[1])
    sample_name = sys.argv[2]
    cpu_num = sys.argv[3]
    lsf_sub(abs_path, sample_name, cpu_num)

