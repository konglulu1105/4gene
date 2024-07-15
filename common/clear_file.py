#!/usr/bin/env python3
import os
import sys
import argparse
import re
import logging
import json
import shutil
BASEPATH = os.path.dirname(os.path.abspath(__file__))

logging.basicConfig(level=logging.INFO,
                    format="[%(asctime)s - %(levelname)s] %(message)s",
                    datefmt="%Y-%m-%d %H:%M:%S")


def parse_arguments():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("--file-list", dest="file_list", nargs="+",
                        help="files that need to be remove")
    return parser.parse_args()


def main():
    args = parse_arguments()
    file_list = args.file_list
    remove_multifiles(file_list)
    
def remove_multifiles(file_list):
    '''
    Remove multiple files

    Args:
    file_list: list of files that need to be remove

    '''
    for file_path in file_list:
        if os.path.isdir(file_path):
            shutil.rmtree(file_path, ignore_errors=True)
        else:
            os.remove(file_path)
        
if __name__ == '__main__':
    main()
