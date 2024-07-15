#!/usr/bin/env python3
import os
import re
import subprocess
import argparse
import logging
from yaml import Loader, Dumper, load
from common.sample_info import SampleInfo, SampleList
from common.base_pipeline import BasePipeline, run_via_shellscript, run_via_snakemake
from common import tools
from collections import OrderedDict

logging.basicConfig(level=logging.INFO,
                    format="[%(asctime)s - %(levelname)s] %(message)s",
                    datefmt="%Y-%m-%d %H:%M:%S")

CONFIG_DICT = dict()
PIPELINE_ROOT = os.path.dirname(os.path.realpath(__file__))
DEFAULT_CONFIG_FILE = f'{PIPELINE_ROOT}/config/four_gene_zhuce.yaml'


def read_config_file(config_file):
    with open(config_file, 'r') as yamlfile:
        return load(yamlfile.read(), Loader=Loader)


def parse_arguments():
    parser = argparse.ArgumentParser(description="main script for analyazing DNA data")
    parser.add_argument("--config-file", "--config", "-c",
                        dest="config_file", default=DEFAULT_CONFIG_FILE, help="config file")
    parser.add_argument("--out-dir", "--path", "-p",
                        required=True,  dest="out_dir")
    parser.add_argument("--via-snakemake", action="store_true", help="if set, will run workflow via snakemake")
    parser.add_argument("--run-mode", required=False, dest="run_mode",
                        choices=["local", "lsf", "none"], default="local",
                        help="run mode, local, lsf or none; \
                                choose local, run in local directly; \
                                choose lsf, submit shell script to lsf cluster; \
                                choose none, will not run the script generated")

    sample_args = parser.add_mutually_exclusive_group(required=True)
    sample_args.add_argument("--sample-name", "--sample", "-s",
                        dest="sample_name", help="sample name")
    sample_args.add_argument("--sample-list", dest="sample_list", help="sample list file")
    parser.add_argument("--read1", "-1", dest="read1", help="fastq read1")
    parser.add_argument("--read2", "-2", dest="read2", help="fastq read2")
    args = parser.parse_args()
    if not args.via_snakemake and args.sample_list is not None:
        raise argparse.ArgumentError(None, "you can only run with sample list file via snakemake by setting via_snakemake to true")
    if args.via_snakemake and args.sample_list is None:
        raise argparse.ArgumentError(None, "when run via snakemake, sample_list must be set")
    if not args.via_snakemake:
        if args.read1 is None or args.read2 is None or args.sample_name is None:
            raise argparse.ArgumentError(None, "you need set fastq_read1, fastq_read2 and sample_name")
    return args


class DNAPipeline(BasePipeline):
    pipeline_name = "dna_pipeline"
    main_workflow = "workflows/DNA/all.smk"
    module_dict = OrderedDict([
        ("trim", "workflows/DNA/trim.smk"),
        ("bwa", "workflows/DNA/bwa.smk"),
        ("varscan", "workflows/DNA/call_varscan.smk"),
        ("anno", "workflows/DNA/anno.smk"),
        ("stat", "workflows/DNA/qc_stat.smk"),
        ("clear", "workflows/DNA/clear.smk")
    ])


def main():
    args = parse_arguments()
    config_file = os.path.abspath(args.config_file)
    global CONFIG_DICT
    CONFIG_DICT = read_config_file(config_file)
    run_mode = args.run_mode
    out_dir = os.path.abspath(args.out_dir)
    if not args.via_snakemake:
        run_via_shellscript(args.sample_name, args.read1, args.read2, out_dir, config_file, run_mode, CONFIG_DICT["snakemake"], PIPELINE_ROOT, DNAPipeline)
    else:
        run_via_snakemake(args.sample_list, out_dir, config_file, run_mode, CONFIG_DICT['snakemake'], PIPELINE_ROOT, DNAPipeline)


if __name__ == '__main__':
    main()
