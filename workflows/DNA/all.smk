import os
import pandas as pd
sample_name_list = pd.read_csv(config["sample_list_file"], sep = '\t', header = 0)['sample_name']
out_dir = config["out_dir"]
base_path = config["base_path"]
python3 = config["python3"]

rule all:
    input:
        expand(os.path.join(out_dir, "TMP", "{sample_name}_step1.done"), sample_name=sample_name_list),
        expand(os.path.join(out_dir, "TMP", "{sample_name}_step2.done"),sample_name=sample_name_list),
        expand(os.path.join(out_dir, "TMP", "{sample_name}_step3.done"),sample_name=sample_name_list),
        expand(os.path.join(out_dir, "TMP", "{sample_name}_step4.done"), sample_name=sample_name_list),
        expand(os.path.join(out_dir, "results", "{sample_name}_SNV.xls"), sample_name=sample_name_list),
        expand(os.path.join(out_dir, "results", "{sample_name}_DNA_QC.xls"), sample_name=sample_name_list)


include: "trim.smk"
include: "bwa.smk"
include: "call_varscan.smk"
include: "anno.smk"
include: "qc_stat.smk"
include: "clear.smk"

