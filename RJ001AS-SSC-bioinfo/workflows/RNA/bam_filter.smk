import os
out_dir = config['out_dir']
python3 = config['python3']
base_path = config['base_path']

rule get_sam1:
    input:
        os.path.join(out_dir, "QC_Bam", "Bam", "{sample_name}_sort.bam")
    output:
        os.path.join(out_dir, "QC_Bam", "Bam", "{sample_name}_sort.sam")
    params:
        samtools = config["samtools"]
    shell:
        """
        {params.samtools} view -h -F 780 \
            {input} \
            -o {output}
        """

rule sam_filter:
    input:
        os.path.join(out_dir, "QC_Bam", "Bam", "{sample_name}_sort.sam")
    output:
        filter_sam = os.path.join(out_dir, "QC_Bam", "Bam", "{sample_name}_sort.sam.filter"),
        reads_stat = os.path.join(out_dir, "TMP", "{sample_name}_reads.stats")
    params:
        python3 = config["python3"],
        base_path = config["base_path"],
        bed = config['RNA']["bed"]
    shell:
        """
        {params.python3} {params.base_path}/RNA/sam_filter.py \
		    --fusion-bed {params.bed} \
			--input-sam {input} \
            --filter-sam {output.filter_sam} \
			--out-stat {output.reads_stat}
        """

rule sam2bam:
    input:
        os.path.join(out_dir, "QC_Bam", "Bam", "{sample_name}_sort.sam.filter")
    output:
        os.path.join(out_dir, "QC_Bam", "Bam", "{sample_name}_sort_filter.bam")
    params:
        samtools = config["samtools"]
    shell:
        """
        {params.samtools} view -b -F 780 \
            -o {output} {input}
        """

rule index_bam:
    input:
        os.path.join(out_dir, "QC_Bam", "Bam", "{sample_name}_sort_filter.bam")
    output:
        os.path.join(out_dir, "QC_Bam", "Bam", "{sample_name}_sort_filter.bam.bai")
    params:
        samtools = config["samtools"]
    shell:
        """
        {params.samtools} index {input}
        """

rule check_step3:
    input:
        os.path.join(out_dir, "results", "{sample_name}_Fusion.xls"),
        os.path.join(out_dir, "QC_Bam", "Bam", "{sample_name}_sort.sam"),
        os.path.join(out_dir, "TMP", '{sample_name}_reads.stats'),
        os.path.join(out_dir, "QC_Bam","Bam",  '{sample_name}_sort.sam.filter')
    output:
        os.path.join(out_dir, "TMP", "{sample_name}_step3.done")
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
            --step-num 3 \
            --process-log {params.out_dir}/TMP/{wildcards.sample_name}_Process.log && \
            touch {params.out_dir}/TMP/{wildcards.sample_name}_step3.done
        """

