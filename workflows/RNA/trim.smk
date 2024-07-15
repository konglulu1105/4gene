import os
out_dir = config["out_dir"]
base_path = config["base_path"]
trimmomatic = config["trimmomatic"]

rule fastqc:
    input:
        r1 = os.path.join(out_dir, "fastq", "{sample_name}.R1.fastq.gz"),
        r2 = os.path.join(out_dir, "fastq", "{sample_name}.R2.fastq.gz")
    output:
        log = os.path.join(out_dir, "TMP", "{sample_name}_fastqc.log")
    params:
        fastqc =  config["fastqc"],
        java = config["java"],
        out_dir = out_dir,
    threads: 4
    shell:
        """
        {params.fastqc} -o {params.out_dir}/TMP \
            -t 4 \
            -java {params.java} \
            {input.r1} \
            {input.r2} \
            2> {output}
        """

rule trimmomatic:
    input:
        r1 = os.path.join(out_dir, "fastq", "{sample_name}.R1.fastq.gz"),
        r2 = os.path.join(out_dir, "fastq", "{sample_name}.R2.fastq.gz")
    output:
        r1 = os.path.join(out_dir, "QC_Bam", "Clean", "{sample_name}_clean_R1.fq.gz"),
        unpaired_r1 = os.path.join(out_dir, "TMP", "{sample_name}_unpaired_R1.fq.gz"),
        r2 = os.path.join(out_dir, "QC_Bam", "Clean", "{sample_name}_clean_R2.fq.gz"),
        unpaired_r2 = os.path.join(out_dir, "TMP", "{sample_name}_unpaired_R2.fq.gz"),
        trim_log = os.path.join(out_dir, "TMP", "{sample_name}_trim.log")
    params:
        java = config["java"],
        trimmomatic = trimmomatic,
        adapters_path = config["DNA"]["adapters_path"],
        out_dir = out_dir
    threads: 2
    shell:
        """
        {params.java} -jar {params.trimmomatic} PE \
            -phred33 \
            -threads 4 \
            {input.r1} \
            {input.r2} \
            {output.r1} \
            {output.unpaired_r1} \
            {output.r2} \
            {output.unpaired_r2} \
            ILLUMINACLIP:{params.adapters_path}:2:30:10:1:true LEADING:15 TRAILING:15 SLIDINGWINDOW:30:25 MINLEN:50 \
            2> {params.out_dir}/TMP/{wildcards.sample_name}_trim.log
        """

rule reseqtools_raw:
    input:
        r1 = os.path.join(out_dir, "fastq", "{sample_name}.R1.fastq.gz"),
        r2 = os.path.join(out_dir, "fastq", "{sample_name}.R2.fastq.gz")
    output:
        os.path.join(out_dir, "TMP", "{sample_name}_fqstat.txt")
    params:
        reseqtools = config["reseqtools"]
    threads: 2
    shell:
        """
        {params.reseqtools} Fqtools stat \
            -InFq {input.r1} \
            -InFq {input.r2} \
            -OutStat {output}
        """

rule reseqtools_clean:
    input:
        r1 = os.path.join(out_dir, "QC_Bam", "Clean", "{sample_name}_clean_R1.fq.gz"),
        r2 = os.path.join(out_dir, "QC_Bam", "Clean", "{sample_name}_clean_R2.fq.gz")
    output:
        os.path.join(out_dir, "TMP", "{sample_name}_clean_fqstat.txt")
    params:
        reseqtools = config["reseqtools"]
    threads: 2
    shell:
        """
        {params.reseqtools} Fqtools stat \
            -InFq {input.r1} \
            -InFq {input.r2} \
            -OutStat {output}
        """
    
rule check_step1:
    input:
        os.path.join(out_dir, "TMP", "{sample_name}_trim.log"),
        os.path.join(out_dir, "TMP", "{sample_name}_fastqc.log"),
        os.path.join(out_dir, "TMP", "{sample_name}_clean_fqstat.txt"),
        os.path.join(out_dir, "TMP", "{sample_name}_fqstat.txt")
    output:
        os.path.join(out_dir, "TMP", "{sample_name}_step1.done")
    params:
        python3 = config["python3"],
        out_dir = out_dir,
        base_path = base_path
    shell:
        """
        {params.python3} {params.base_path}/common/check_step.py \
            --scan-path {params.out_dir} \
            --sample-name {wildcards.sample_name} \
            --workflow-name RNA \
            --step-num 1 \
            --process-log {params.out_dir}/TMP/{wildcards.sample_name}_Process.log && \
        touch {params.out_dir}/TMP/{wildcards.sample_name}_step1.done
        """
