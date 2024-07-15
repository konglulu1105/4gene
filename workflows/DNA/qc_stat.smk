import os

out_dir = config["out_dir"]
python3 = config["python3"]
base_path = config["base_path"]

# cmd1
rule DepthofCoverage:
    input:
        sort_bam = os.path.join(out_dir, "QC_Bam", "Bam", "{sample_name}_sort.bam"),
        bam_index = os.path.join(out_dir, "QC_Bam", "Bam", "{sample_name}_sort.bam.bai")
    output:
        depth_summary_file = os.path.join(out_dir, "TMP", "{sample_name}_sort_DepthOfCoverage.sample_summary"),
        depth_base_file = os.path.join(out_dir, "TMP", "{sample_name}_sort_DepthOfCoverage")
    params:
        java = config["java"],
        gatk3 = config["GATK3"],
        out_dir = out_dir,
        ref = config["DNA"]["ref"],
        bed = config["DNA"]["bed"]
    threads: 4
    shell:
        """
        {params.java} -jar {params.gatk3} \
		    -T DepthOfCoverage \
	        --start 50 \
		    --stop 300000 \
			--filter_mismatching_base_and_quals \
			-R {params.ref} \
			-o {params.out_dir}/TMP/{wildcards.sample_name}_sort_DepthOfCoverage \
			-I {params.out_dir}/QC_Bam/Bam/{wildcards.sample_name}_sort.bam \
			-L {params.bed} \
			-ct 100 -ct 200 -ct 500
        """
# cmd3
rule CollectInsertSizeMetrics:
    input:
        sort_bam = os.path.join(out_dir, "QC_Bam", "Bam", "{sample_name}_sort.bam"),
        bam_index = os.path.join(out_dir, "QC_Bam", "Bam", "{sample_name}_sort.bam.bai")
    output:
        insert_size_metrics = os.path.join(out_dir, "TMP", "{sample_name}_insert_size_metrics.txt")
    params:
        java = config["java"],
        picard = config["picard"],
        out_dir = out_dir
    threads: 2
    shell:
        """
        {params.java} -jar {params.picard} CollectInsertSizeMetrics \
            I={params.out_dir}/QC_Bam/Bam/{wildcards.sample_name}_sort.bam \
			O={params.out_dir}/TMP/{wildcards.sample_name}_insert_size_metrics.txt \
			H={params.out_dir}/TMP/{wildcards.sample_name}_insert_size_histogram.pdf \
			M=0.5
        """

# cmd4
rule bedtools_intersect:
    input:
        sort_bam = os.path.join(out_dir, "QC_Bam", "Bam", "{sample_name}_sort.bam"),
        bam_index = os.path.join(out_dir, "QC_Bam", "Bam", "{sample_name}_sort.bam.bai")
    output:
        sort_cover_bam = os.path.join(out_dir, "QC_Bam", "Bam", "{sample_name}_sort_cover.bam")
    params:
        bedtools = config["bedtools"],
        out_dir = out_dir,
        bed = config["DNA"]["bed"]
    threads: 2
    shell:
        """
        {params.bedtools} intersect \
            -wa \
            -a {params.out_dir}/QC_Bam/Bam/{wildcards.sample_name}_sort.bam \
			-b {params.bed} \
			> {params.out_dir}/QC_Bam/Bam/{wildcards.sample_name}_sort_cover.bam
        """
# cmd6
rule sort_bam_flagstat:
    input:
        sort_bam = os.path.join(out_dir, "QC_Bam", "Bam", "{sample_name}_sort.bam"),
        bam_index = os.path.join(out_dir, "QC_Bam", "Bam", "{sample_name}_sort.bam.bai")
    output:
        flag_stat_file = os.path.join(out_dir, "TMP", "{sample_name}_sort_flagstat.txt")
    params:
        samtools = config["samtools"],
        out_dir = out_dir
    threads: 2
    shell:
        """
        {params.samtools} flagstat \
            {params.out_dir}/QC_Bam/Bam/{wildcards.sample_name}_sort.bam \
            > {params.out_dir}/TMP/{wildcards.sample_name}_sort_flagstat.txt
        """
# cmd7
rule sort_cover_bam_flagstat:
    input:
        sort_cover_bam = os.path.join(out_dir, "QC_Bam", "Bam", "{sample_name}_sort_cover.bam"),
        bam_index = os.path.join(out_dir, "QC_Bam", "Bam", "{sample_name}_sort.bam.bai")
    output:
        flag_stat_file = os.path.join(out_dir, "TMP", "{sample_name}_sort_cover_flagstat.txt")
    params:
        samtools = config["samtools"],
        out_dir = out_dir
    threads: 2
    shell:
        """
        {params.samtools} flagstat \
            {params.out_dir}/QC_Bam/Bam/{wildcards.sample_name}_sort_cover.bam \
            > {params.out_dir}/TMP/{wildcards.sample_name}_sort_cover_flagstat.txt
        """

# cmd8
rule gather_qc_stat:
    input:
        raw_fqstat_file = os.path.join(out_dir, "TMP", "{sample_name}_fqstat.txt"),
        clean_fqstat_file = os.path.join(out_dir, "TMP", "{sample_name}_clean_fqstat.txt"),
        sort_flagstat = os.path.join(out_dir, "TMP", "{sample_name}_sort_flagstat.txt"),
        sort_cover_flagstat = os.path.join(out_dir, "TMP", "{sample_name}_sort_cover_flagstat.txt"),
        insert_metrics = os.path.join(out_dir, "TMP", "{sample_name}_insert_size_metrics.txt"),
        depth_base_file =  os.path.join(out_dir, "TMP", "{sample_name}_sort_DepthOfCoverage"),
        depth_summary_file = os.path.join(out_dir, "TMP", "{sample_name}_sort_DepthOfCoverage.sample_summary")
    output:
        os.path.join(out_dir, "results", "{sample_name}_DNA_QC.xls")
    params:
        python3 = config["python3"],
        base_path = config["base_path"]
    shell:
        """
        {params.python3} {params.base_path}/DNA/dna_qc.py \
            --sample_name {wildcards.sample_name} \
            --raw_fqstat_file {input.raw_fqstat_file} \
            --clean_fqstat_file {input.clean_fqstat_file} \
            --sort_flagstat {input.sort_flagstat} \
            --sort_cover_flagstat {input.sort_cover_flagstat} \
            --insert_metrics {input.insert_metrics} \
            --depth_base_file {input.depth_base_file} \
            --depth_summary_file {input.depth_summary_file} \
            --out_file {output}
        """

# cmd 9
rule check_step4:
    input:
        os.path.join(out_dir, "results", "{sample_name}_DNA_QC.xls")
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
            --workflow-name DNA \
            --step-num 4 \
            --process-log {params.out_dir}/TMP/{wildcards.sample_name}_Process.log && \
        touch {params.out_dir}/TMP/{wildcards.sample_name}_step4.done
        """
