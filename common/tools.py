import sys
import os
import re
import subprocess
import logging
logging.basicConfig(level=logging.INFO,
                    format="[%(asctime)s - %(levelname)s] %(message)s",
                    datefmt="%Y-%m-%d %H:%M:%S")

def run_cmd(cmd_args):
    """run input command through subprocess.check_call

    Args:
        cmd_args (string): command args that need to be executed
    """
    logging.info(f"CMD: {cmd_args}")
    subprocess.check_call(cmd_args, shell=True)
    
def run_shell_script(shell_script):
    """run shell script thougn subprocess.check_call

    Args:
        shell_script (string): path to the shell script that need to run
    """
    logging.info(f"shell script: {shell_script}")
    cmd_args = " ".join(["bash", shell_script])
    subprocess.check_call(cmd_args, shell=True)
    
def create_dir(dir_name):
    """create directory with given name, if not exists

    Args:
        dir_name (string): path to directory that need to be created

    Returns:
        string: path to the given name of directory
    """
    logging.info(f"creating {dir_name}")
    if not os.path.exists(dir_name):
        try:
            os.makedirs(dir_name, 0o755)
        except FileExistsError:
            logging.info(f"{dir_name} already exists")
    return dir_name


def tsv2list(tsv_file):
    """to load tsv file into list of dict (key=header from 1st line, value )

    Args:
        tsv_file (string): path to tsv file
        
    Return:
        list of dict
    
    Raises:
        raise AssertionError if a line has different number of columns from the header line
    """
    result_list = []
    headers = list()
    with open(tsv_file, 'r') as finp:
        headers = finp.readline().strip("\n").split("\t")
        idx = 1
        for line in finp:
            idx += 1
            cols = line.strip("\n").split("\t")
            assert len(headers) == len(cols), \
                f"{tsv_file} L{idx} has different number of columns {len(cols)} but header line {len(headers)}"
            result_list.append(dict(zip(headers, cols)))
    return headers, result_list
