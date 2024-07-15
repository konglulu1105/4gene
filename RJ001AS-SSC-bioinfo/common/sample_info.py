#!/usr/bin/env python3
import os
from typing import List
try:
    import tools
except ImportError:
    from . import tools
class SampleInfo(object):
    def __init__(self, sample_name: str,  out_dir: str, fastq_r1: str="", fastq_r2: str="") -> None:
        self.sample_name = sample_name
        self.fastq_r1 = fastq_r1
        self.fastq_r2 = fastq_r2
        self.out_dir = out_dir
        self.run_script = os.path.join(out_dir, "Run_scripts", f"{self.sample_name}.sh")
        self.run_script_dir = os.path.join(out_dir, "Run_scripts")

    def create_out_folders(self):
        tools.create_dir(f"{self.out_dir}/TMP")
        tools.create_dir(f"{self.out_dir}/Varscan")
        tools.create_dir(f"{self.out_dir}/Varscan/VCF")
        tools.create_dir(f"{self.out_dir}/results")
        tools.create_dir(f"{self.out_dir}/Run_scripts")
        tools.create_dir(f"{self.out_dir}/QC_Bam")
        tools.create_dir(f"{self.out_dir}/QC_Bam/Bam")
        tools.create_dir(f"{self.out_dir}/QC_Bam/Clean")
        tools.create_dir(f"{self.out_dir}/fastq")

    def get_link_fastq(self):
        fastq_dir = f"{self.out_dir}/fastq"
        os.system(f"ln -sfv {self.fastq_r1} {fastq_dir}/{self.sample_name}.R1.fastq.gz")
        os.system(f"ln -sfv {self.fastq_r2} {fastq_dir}/{self.sample_name}.R2.fastq.gz")
        
    def get_module_script_path(self, module: str):
        return f"{self.run_script_dir}/{self.sample_name}_{module}.sh"
    
    
class SampleList(object):
    def __init__(self, sample_list_file: str, out_dir: str) -> None:
        self.sample_list_file = sample_list_file
        self.out_dir = out_dir
        self.samples = self._read_sample_list()
        for sample_info in self.samples:
            sample_info.create_out_folders()
            sample_info.get_link_fastq()
        
    def _read_sample_list(self) -> List[SampleInfo]:
        headers, result_list = tools.tsv2list(self.sample_list_file)
        samples = list()
        for result in result_list:
            sample_info = SampleInfo(
                sample_name=result["sample_name"],
                out_dir=self.out_dir,
                fastq_r1=result["fastq_r1"],
                fastq_r2=result["fastq_r2"]
            )
            samples.append(sample_info)
        return samples
            
        