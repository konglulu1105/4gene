#!/usr/bin/env python3
import os
import sys
import argparse
import re
import logging
import json
BASEPATH = os.path.dirname(os.path.abspath(__file__))

logging.basicConfig(level=logging.INFO,
                    format="[%(asctime)s - %(levelname)s] %(message)s",
                    datefmt="%Y-%m-%d %H:%M:%S")


def parser_arguments():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("--scan-path", dest="scan_path",
                        help="scan path, result files lies")
    parser.add_argument(
        "--sample-name", dest="sample_name", help="sample name")
    parser.add_argument("--workflow-name", dest="workflow_name", choices=[
                        "common", "DNA", "RNA"], default="common", help="name of workflow")
    parser.add_argument("--step-num", dest="step_num",
                        required=True, help="num of step")
    parser.add_argument("--process-log", dest="process_log",
                        required=True, help="process log")
    return parser.parse_args()


def main():
    '''
    Select checker based on the workflow_name,check the analysis step,and record to log

    '''
    args = parser_arguments()
    scan_path = args.scan_path
    sample_name = args.sample_name
    workflow_name = args.workflow_name
    step_num = args.step_num
    process_log = args.process_log
    checker_class = None
    if workflow_name == "DNA":
        checker_class = DNAStepchecker
    elif workflow_name == "RNA":
        checker_class = RNAStepchecker
    else:
        checker_class = WfStepChecker
    step_checker = checker_class(scan_path=scan_path,
                                 sample_name=sample_name,
                                 process_log=process_log)
    if step_num in checker_class.steps:
        step_checker.check_step(step_num)
    else:
        raise ValueError(f"step checker cannot check this {step_num}")


class WfStepChecker(object):
    steps = {'1': ['_fastqc.log', '_trim.log', '_fqstat.txt', '_clean_fqstat.txt'],
             '2': ['_bwa.log', '_sort.bam', '_sort.bam.bai']}

    def __init__(self, scan_path, sample_name, process_log) -> None:
        self.scan_path = scan_path
        self.sample_name = sample_name
        self.process_log = process_log

    def find_file(self):
        '''
        Find the path of result base on sample name

        Returns:
            the path of result
        '''
        result = []
        for root, list, files in os.walk(self.scan_path):
            for file in files:
                if self.sample_name in file:
                    write = os.path.join(root, file)
                    result.append(write)
        return result

    def check_step(self, step_num):
        '''
        Check if the analysis is complete

         Args:
            step_num:the num of analysis steps
        '''
        sample_results = self.find_file()
        step = self.steps[step_num]
        total = []
        for i in range(len(step)):
            for j in range(len(sample_results)):
                if sample_results[j].endswith(step[i]) and os.path.getsize(sample_results[j]) != 0:
                    total.append(i)
        List = list(set(total))
        if len(List) == len(step):
            flag = '1'
        else:
            flag = '0'
        self.update_process_log(step_num, flag)

    def update_process_log(self, step_num, flag):
        with open(self.process_log, 'a') as f:
            f.write("%s\n" % ('step' + step_num + ':' + flag))


class DNAStepchecker(WfStepChecker):
    steps = {'1': ['_fastqc.log', '_trim.log', '_fqstat.txt', '_clean_fqstat.txt'],
             '2': ['_bwa.log', '_sort.bam', '_sort.bam.bai'],
             '3': ['_sort.mpileup', '_varscan_snp.vcf', '_varscan_indel.vcf', '_combined.vcf'],
             '4': ['_SNV.xls', '_sort_DepthOfCoverage', '_insert_size_metrics.txt', '_sort_cover.bam',
                   '_sort_flagstat.txt', '_sort_cover_flagstat.txt', '_DNA_QC.xls']}

    def update_process_log(self, step_num, flag):
        with open(self.process_log, 'a') as f:
            f.write("%s\n" % ('DNA step' + step_num + ':' + flag))


class RNAStepchecker(WfStepChecker):
    steps = {'1': ['_fastqc.log', '_trim.log', '_fqstat.txt', '_clean_fqstat.txt'],
             '2': ['_bwa.log', '_sort.bam', '_sort.bam.bai'],
             '3': ['_sort.sam', '_reads.stats', '_sort.sam.filter', '_sort_filter.bam'],
             '4': ['_Fusion.xls', '_sort_flagstat.txt', '_RNA_QC.xls']}

    def update_process_log(self, step_num, flag):
        with open(self.process_log, 'a') as f:
            f.write("%s\n" % ('RNA step' + step_num + ':' + flag))


if __name__ == '__main__':
    main()
