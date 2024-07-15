#!/usr/bin/env python3
import sys
import os
import time
import re
import subprocess
import logging
from snakemake import io, utils
from snakemake import workflow as _workflow

class SMKRunner(object):
    def __init__(self, snakemake_bin, base_path, snakefile, out_dir, config_file) -> None:
        self.snakemake_bin = snakemake_bin
        self.base_path = base_path
        self.snakefile = snakefile
        self.out_dir = out_dir
        self.config_file = config_file
        
    def get_cmd_args(self, sample_list_file: str):
        cmd_args = []
        cmd_args.extend([self.snakemake_bin])
        cmd_args.extend(["-s", self.snakefile])
        cmd_args.extend(["--directory", self.out_dir])
        cmd_args.extend(["--configfile", self.config_file])
        cmd_args.extend(["--config", f"sample_list_file={sample_list_file}", f"out_dir={self.out_dir}", f"base_path={self.base_path}"])
        cmd_args.extend(["--cores"])
        cmd_args.extend(["2>&1"])
        cmd_args.extend(["|", "tee", "-a", f"{self.out_dir}/snakemake.log",])
        return cmd_args
        
    def get_shellcmd_list(self, sample_name):
        config_dict = {
            "sample_name": sample_name,
            "out_dir": self.out_dir,
            "base_path": self.base_path
        }
        overwrite_config = dict()
        utils.update_config(overwrite_config, config_dict)
        utils.update_config(overwrite_config, io.load_configfile(self.config_file))
        workflow = _workflow.Workflow(self.snakefile, overwrite_config=overwrite_config)
        workflow.include(self.snakefile)
        workflow.check()
        wildcards = io.Wildcards()
        wildcards.append(sample_name)
        wildcards._add_name("sample_name")
        shellcmd_list = list()
        for rule in workflow.rules:
            if not rule.is_shell:
                continue
            if rule.name.startswith("disable_"):
                continue
            raw_string = rule.shellcmd.format(input=rule.input, output=rule.output, 
                                              params=rule.params, sample_name=sample_name, 
                                              wildcards=wildcards)
            raw_string = re.sub('\s+',' ', raw_string).strip()
            shellcmd_list.append(raw_string.format(sample_name=sample_name).lstrip().strip("\n"))
        return shellcmd_list
