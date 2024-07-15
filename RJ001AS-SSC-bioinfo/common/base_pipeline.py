#!/usr/bin/env python3
import sys
import os
from collections import OrderedDict
from snakemake import io, utils
import logging
try:
    from .base_runner import SMKRunner
    from .sample_info import SampleInfo, SampleList
    from . import tools
except ImportError:
    from base_runner import SMKRunner
    from sample_info import SampleInfo, SampleList
    import tools
logging.basicConfig(level=logging.INFO,
                    format="[%(asctime)s - %(levelname)s] %(message)s",
                    datefmt="%Y-%m-%d %H:%M:%S")

class BasePipeline(object):
    pipeline_name = "" 
    module_dict =  OrderedDict()
    main_workflow = ""
    
    def __init__(self, base_path: str, config_file: str, snakemake_bin: str ="") -> None:
        self.snakemake_bin = snakemake_bin
        self.base_path = base_path
        self.config_file = config_file
        self.config = io.load_configfile(self.config_file)
    
    def get_snakemake_run_args(self, sample_list: SampleList):
        smk = os.path.join(self.base_path, self.main_workflow)
        smk_runner = SMKRunner(self.snakemake_bin, self.base_path, smk, sample_list.out_dir, self.config_file)
        return smk_runner.get_cmd_args(sample_list.sample_list_file)

    def generate_shell_script(self, sample_info: SampleInfo):
        sample_run_script = sample_info.run_script
        outf = open(sample_info.run_script, 'w')
        for module_name, smk in self.module_dict.items():
            shellcmd_list = self.generate_module_script(sample_info, module_name, smk)
            with open(sample_info.get_module_script_path(module_name), "w") as fout:
                fout.write("\n".join(shellcmd_list) + "\n")
            outf.write(f"\n#################module {module_name}##################\n")
            outf.write("\n".join(shellcmd_list) + "\n")
        outf.close()
        return sample_run_script

    def generate_module_script(self, sample_info: SampleInfo, module_name: str, smk: str):
        smk = os.path.join(self.base_path, smk)
        smk_runner = SMKRunner(self.snakemake_bin, self.base_path, smk, sample_info.out_dir, self.config_file)
        try:
            shellcmd_list = smk_runner.get_shellcmd_list(sample_info.sample_name)
        except Exception as e:
            logging.error(f"ERROR FOUND in module {module_name}: {smk}, " + \
                          "when generating shell script for {sample_info.sample_name}")
            raise Exception
        return shellcmd_list
    
    def process_run_script(self, sample_info: SampleInfo, run_script: str, run_mode: str):
        if run_mode == "lsf":
            self.submit_script_to_lsf(sample_info, run_script)
        elif run_mode == "local":
            tools.run_shell_script(run_script)
        else:
            logging.info("run mode choose none, will not run the script generated")
            pass
        
    def submit_script_to_lsf(self, sample_info: SampleInfo, sample_run_script: str):
        sample_name = sample_info.sample_name
        shell_script = sample_run_script
        print(shell_script)
        if not os.path.exists(shell_script):
            raise Exception(f"cannot find {shell_script}, please check it")
        submit_args = " ".join([self.config["python3"], f"{self.base_path}/common/lsf.py",
                                shell_script, sample_name, "4"])
        print(submit_args)
        logging.info(f"submit shell script {shell_script} to lsf")
        tools.run_cmd(submit_args)
        
def run_via_shellscript(sample_name, fastq_read1, fastq_read2, out_dir, config_file, run_mode, snakemake_bin, base_path, pipe_class):
    fastq_read1 = os.path.abspath(fastq_read1)
    fastq_read2 = os.path.abspath(fastq_read2)
    sample_info = SampleInfo(sample_name=sample_name,
                             fastq_r1=fastq_read1,
                             fastq_r2=fastq_read2,
                             out_dir=out_dir)
    sample_info.create_out_folders()
    sample_info.get_link_fastq()
    pipe = pipe_class(base_path, config_file=config_file, snakemake_bin=snakemake_bin)
    sample_run_script = pipe.generate_shell_script(sample_info)
    pipe.process_run_script(sample_info, sample_run_script, run_mode)

def run_via_snakemake(sample_list_file, out_dir, config_file, run_mode, snakemake_bin, base_path, pipe_class):
    sample_list_file = os.path.abspath(sample_list_file)
    pipe = pipe_class(base_path, config_file=config_file, snakemake_bin=snakemake_bin)
    sample_list = SampleList(sample_list_file, out_dir)
    cmd_args = pipe.get_snakemake_run_args(sample_list)
    tools.run_cmd(" ".join(cmd_args))