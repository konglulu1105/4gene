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
    parser.add_argument("--fusion-raw", dest="fusion_raw", required=True,
                        help="fusion raw")
    parser.add_argument("--data-qc", dest="data_qc", required=True,
                        help="qc data for RNA seq")
    parser.add_argument("--out-file", dest="out_file",
                        required=True, help="out file")
    return parser.parse_args()


def main():
    args = parse_arguments()
    fusion_raw = args.fusion_raw
    data_qc = args.data_qc
    out_file = args.out_file
    filter_fusion(fusion_raw, data_qc, out_file)

def load_data_qc(data_qc_file):
    with open(data_qc_file) as finp:
        headers = finp.readline().strip("\n").split("\t")
        line = finp.readline().strip("\n").split("\t")
        return dict(zip(headers, line))


def filter_fusion(fusion_raw, data_qc_file, out_file):
    """filter fusion result with internal rules

    Args:
        fusion_raw (string): path to the raw fusion result file
        data_qc_file (string): path to the RNA qc data file
        out_file (string): path to output file, contained fusion after filtering
    """
    clean_data = load_data_qc(data_qc_file)
    fout = open(out_file, 'w')
    with open(fusion_raw) as finp:
        title = next(finp)
        fout.write(title)
        for line in finp:
            line = line.strip("\n")
            line = line.split("\t")
            support = int(line[6])
            if clean_data == 0:
                norm = support
            else:
                norm = int(support * 50 / float(clean_data['Clean_Data(Mb)']))
            if norm >= 0:
                line[6] = str(norm)
                line[7] = "Pass"
                fout.write("\t".join(line) + "\n")
    fout.close()


if __name__ == '__main__':
    main()