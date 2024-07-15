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


def parse_arguments():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("--sites-json", dest="sites_json",
                        help="sites_json, if not specified, use common/sites.json by default",
                        default=f"{BASEPATH}/sites.json")
    parser.add_argument("--inp-file", dest="inp_file",
                        required=True, help="path to the input file")
    parser.add_argument("--out-file", dest="out_file",
                        required=True, help="path to the output file")
    parser.add_argument("--sample-name", dest="sample_name", required=True, help="sample name")
    parser.add_argument("--variant-type", dest="variant_type",
                        required=True, choices=["SNV", "FUSION"], help="variant type")
    return parser.parse_args()

def main():
    args = parse_arguments()
    sites_json = args.sites_json
    inp_file = args.inp_file
    out_file = args.out_file
    variant_type = args.variant_type
    sample_name = args.sample_name
    result_processor = ResultProcessor(sites_json)
    result_processor.execute(inp_file, out_file, variant_type, sample_name)

class ResultProcessor(object):
    '''
    Summarize the input file and mutation types into the output file


    '''
    def __init__(self, sites_json, log=logging) -> None:
        self.sites_json = sites_json
        self.log = log
        self.load_sites_json()

    def load_sites_json(self):
        inf = open(self.sites_json, 'r')
        json_dict = json.load(inf)
        self.snv_sites = json_dict["SNV"]
        self.fusions = json_dict["FUSION"]

    def execute(self, inp_file, out_file, variant_type, sample_name):
        '''
        Summarize and output information based on mutation types

         Args:
            inp_file:the path of input file
            out_file:the path of output file
            variant_type:the type of variant
            sample_name:the name of sample
        '''
        headers, result_list = load_inp_file(inp_file)
        if variant_type == "SNV":
            key_colnames = ['Chr', 'Position_Start', 'Reference', 'Mutation']
            sites = self.snv_sites
        elif variant_type == "FUSION":
            key_colnames = ['Gene']
            sites = self.fusions
        else:
            raise Exception(
                f"cannot recogize the variant_type {variant_type} you specified")
        out_result_list = self.add_extra_sites(sample_name,
            headers, result_list, sites, key_colnames)
        fout = open(out_file, 'w')
        fout.write('\t'.join(headers) + '\n')
        for result in out_result_list:
            fout.write('\t'.join([result[h] for h in headers]) + "\n")
        fout.close()

    @staticmethod
    def add_extra_sites(sample_name, headers, result_list, sites, col_names):
        '''
        Append the massage of site to result list

        Args:
            sample_name:the name of sample
            headers:the headers of result list
            result_list:the list with sample massage
            sites:list with different site
            col_names:the column name

        Returns:
            list of all sample result
        '''
        out_result_list = list()
        result_dict = dict()
        for result in result_list:
            key = tuple(result[col_name] for col_name in col_names)
            result_dict[key] = result
        for site in sites:
            key = tuple(site[col_name] for col_name in col_names)
            if key in result_dict:
                out_result_list.append(result_dict[key])
            else:
                out_dict = dict()
                out_dict['Sample_Name'] = sample_name
                out_dict['Status'] = "Fail"
                out_dict['Support_Reads'] = 0
                out_dict.update(site)
                out_result_list.append(
                    {h: str(out_dict.get(h, '/')) for h in headers})
        return out_result_list


def load_inp_file(inp_file):
    '''
    Load input file content

    Args:
        inp_file：absolute path of input file

    Returns:
        list of input file with header and remaining content
    '''
    result_list = []
    headers = []
    with open(inp_file, 'r') as finp:
        headers = finp.readline().strip().strip('\n').split('\t')
        idx = 1
        for line in finp.readlines():
            idx += 1
            cols = line.strip('\n').split('\t')
            assert len(headers) == len(cols), \
                f"{inp_file} L{idx} has different number of columns {len(cols)} but header line {len(headers)}"
            result_list.append(dict(zip(headers, cols)))
    return headers, result_list


if __name__ == '__main__':
    main()