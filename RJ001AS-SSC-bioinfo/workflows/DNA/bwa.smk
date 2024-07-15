import os
out_dir = config["out_dir"]
python3 = config["python3"]
base_path = config["base_path"]

rule bwa:
    input:
        r1 = os.path.join(out_dir, "QC_Bam", "Clean", "{sample_name}_clean_R1.fq.gz"),
        r2 = os.path.join(out_dir, "QC_Bam", "Clean", "{sample_name}_clean_R2.fq.gz")
    output:
        bam =  os.path.join(out_dir, "QC_Bam", "Bam", "{sample_name}.bam"),
        log = os.path.join(out_dir, "TMP", "{sample_name}_bwa.log")
    params:
        samtools = config["samtools"],
        bwa = config["bwa"],
        ref = config["DNA"]["ref"],
        out_dir = out_dir
    threads: 4
    shell:
        """
        {params.bwa} mem -t 4 \
            -R '@RG\\tID:{wildcards.sample_name}\\tSM:{wildcards.sample_name}\\tPL:illumina\\tLB:library\\tPU:PE' \
            {params.ref} {input.r1} {input.r2} \
            > {output.bam} \
            2> {params.out_dir}/TMP/{wildcards.sample_name}_bwa.log
        """

rule sort:
    input:
        bam =  os.path.join(out_dir, "QC_Bam", "Bam", "{sample_name}.bam")
    output:
        sort_bam = os.path.join(out_dir, "QC_Bam", "Bam", "{sample_name}_sort.bam"),
    params:
        samtools = config["samtools"]
    threads: 5
    shell:
        """
        {params.samtools} sort --threads 5 {input} -o {output.sort_bam}
        """

rule bam_index:
    input:
        sort_bam = os.path.join(out_dir, "QC_Bam", "Bam", "{sample_name}_sort.bam")
    output:
        bam_index = os.path.join(out_dir, "QC_Bam", "Bam", "{sample_name}_sort.bam.bai")
    params:
        samtools = config["samtools"]
    threads: 2
    shell:
        """
        {params.samtools} index {input.sort_bam}
        """

rule check_step:
    input:
        sort_bam = os.path.join(out_dir, "QC_Bam", "Bam", "{sample_name}_sort.bam"),
        bam_index = os.path.join(out_dir, "QC_Bam", "Bam", "{sample_name}_sort.bam.bai"),
        bwa_log =  os.path.join(out_dir, "TMP", "{sample_name}_bwa.log")
    output:
        os.path.join(out_dir, "TMP", "{sample_name}_step2.done")
    params:
        python3 = python3,
        out_dir = out_dir,
        base_path = base_path
    shell:
        """
        {params.python3} {params.base_path}/common/check_step.py \
            --scan-path {params.out_dir} \
            --sample-name {wildcards.sample_name} \
            --workflow-name DNA \
            --step-num 2 \
            --process-log {params.out_dir}/TMP/{wildcards.sample_name}_Process.log && \
        touch {params.out_dir}/TMP/{wildcards.sample_name}_step2.done
        """
