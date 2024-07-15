import os
out_dir = config['out_dir']
python3 = config['python3']
base_path = config["base_path"]

rule flagstat:
    input:
        os.path.join(out_dir, "QC_Bam", "Bam", "{sample_name}_sort.bam")
    output:
        os.path.join(out_dir, "TMP", "{sample_name}_sort_flagstat.txt")
    params:
        samtools = config["samtools"]
    threads: 2
    shell:
        """
        {params.samtools} flagstat \
            {input} \
            > {output}
        """

rule get_stat:
    input:
        fqstat_raw = os.path.join(out_dir, "TMP", "{sample_name}_fqstat.txt"),
        fqstat_clean = os.path.join(out_dir, "TMP", "{sample_name}_clean_fqstat.txt"),
        flagstat =  os.path.join(out_dir, "TMP", "{sample_name}_sort_flagstat.txt"),
        reads_stat = os.path.join(out_dir, "TMP", "{sample_name}_reads.stats")
    output:
        os.path.join(out_dir, "results", "{sample_name}_RNA_QC.xls")
    params:
        python3 = config["python3"],
        base_path = config["base_path"],
        bed = config['RNA']["bed"]
    shell:
        """
        {params.python3} {params.base_path}/RNA/rna_qc.py \
            --sample-name {wildcards.sample_name} \
            --fqstat-raw {input.fqstat_raw} \
            --fqstat-clean {input.fqstat_clean} \
            --flagstat {input.flagstat} \
            --reads-stat {input.reads_stat} \
			--out-file {output}
        """


rule check_fusion_qc:
    input:
        fusion_raw = os.path.join(out_dir, "TMP", "{sample_name}_fusion_raw.txt"),
        rna_qc = os.path.join(out_dir, "results", "{sample_name}_RNA_QC.xls"),
    output:
        os.path.join(out_dir, "TMP", "{sample_name}_Fusion_filter.xls")
    params:
        python3 = config["python3"],
        base_path = config['base_path']
    shell:
        """
        {params.python3} {params.base_path}/RNA/fusion_qc.py \
            --fusion-raw {input.fusion_raw} \
            --data-qc {input.rna_qc} \
            --out-file {output}
        """

rule add_extra_fusion:
    input:
        os.path.join(out_dir, "TMP", "{sample_name}_Fusion_filter.xls")
    output:
        os.path.join(out_dir, "results", "{sample_name}_Fusion.xls")
    params:
        python3 = python3,
        base_path = base_path,
        out_dir = out_dir
    shell:
        """
        {params.python3} {params.base_path}/common/add_extra_sites.py \
            --inp-file {input} \
            --out-file {output} \
            --variant-type FUSION \
            --sample-name {wildcards.sample_name}
        """

rule check_step4:
    input:
        os.path.join(out_dir, "results", "{sample_name}_Fusion.xls"),
        os.path.join(out_dir, "TMP", "{sample_name}_sort_flagstat.txt"),
        os.path.join(out_dir, "results", "{sample_name}_RNA_QC.xls")
    output:
        os.path.join(out_dir, "TMP", "{sample_name}_step4.done")
    params:
        python3 = python3,
        out_dir = out_dir,
        base_path = base_path
    shell:
        """
        {params.python3} {params.base_path}/common/check_step.py \
            --scan-path {params.out_dir} \
            --sample-name {wildcards.sample_name} \
            --workflow-name RNA \
            --step-num 4 \
            --process-log {params.out_dir}/TMP/{wildcards.sample_name}_Process.log && \
        touch {params.out_dir}/TMP/{wildcards.sample_name}_step4.done
        """


        




