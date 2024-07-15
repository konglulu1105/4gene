import os
out_dir = config["out_dir"]
base_path = config["base_path"]
python3 = config["python3"]

# cmd 2
rule sort_bam_mpileup:
    input:
        bam = os.path.join(out_dir, "QC_Bam", "Bam", "{sample_name}_sort.bam")
    output:
        mpileup_file = os.path.join(out_dir, "QC_Bam", "Bam", "{sample_name}_sort.mpileup")
    params:
        samtools = config["samtools"],
        out_dir = out_dir,
        ref = config["DNA"]["ref"],
        bed = config["DNA"]["bed"],
        target_bed = config["DNA"]["target_bed"]
    threads: 4
    shell:
        """
        {params.samtools} mpileup \
            -d 0 \
            -A -B \
            -f {params.ref} \
			-l {params.target_bed} \
			{params.out_dir}/QC_Bam/Bam/{wildcards.sample_name}_sort.bam > \
			{params.out_dir}/QC_Bam/Bam/{wildcards.sample_name}_sort.mpileup
        """
# cmd 4
rule mpileup2snp:
    input:
        mpileup_file = os.path.join(out_dir, "QC_Bam", "Bam", "{sample_name}_sort.mpileup")
    output:
        snp_vcf = os.path.join(out_dir, "Varscan", "VCF", "{sample_name}_varscan_snp.vcf")
    params:
        java = config["java"],
        varscan = config["varscan"],
        out_dir = out_dir
    threads: 4
    shell:
        """
        {params.java} -jar {params.varscan} mpileup2snp \
            {params.out_dir}/QC_Bam/Bam/{wildcards.sample_name}_sort.mpileup \
            --min-coverage 8 \
			--min-reads2 1 \
            --min-var-freq 0.0001 \
            --min-freq-for-hom 0.90 \
            --p-value 0.99 \
            --strand-filter 0 \
            --output-vcf 1 > \
            {params.out_dir}/Varscan/VCF/{wildcards.sample_name}_varscan_snp.vcf 2> \
			{params.out_dir}/Varscan/VCF/{wildcards.sample_name}_varscan_snp.err
        """

# cmd 6
rule mpileup2indel:
    input:
        mpileup_file = os.path.join(out_dir, "QC_Bam", "Bam", "{sample_name}_sort.mpileup")
    output:
        indel_vcf = os.path.join(out_dir, "Varscan", "VCF", "{sample_name}_varscan_indel.vcf")
    params:
        java = config["java"],
        varscan = config["varscan"],
        out_dir = out_dir
    threads: 4
    shell:
        """
        {params.java} -jar {params.varscan} mpileup2indel \
            {params.out_dir}/QC_Bam/Bam/{wildcards.sample_name}_sort.mpileup \
            --min-coverage 8 \
			--min-reads2 1 \
			--min-var-freq 0.0001 \
			--min-freq-for-hom 0.90 \
			--p-value 0.99 \
			--strand-filter 0 \
			--output-vcf 1 > \
			{params.out_dir}/Varscan/VCF/{wildcards.sample_name}_varscan_indel.vcf 2> \
			{params.out_dir}/Varscan/VCF/{wildcards.sample_name}_varscan_indel.err
        """

# cmd 11
rule variant_combine_vcf:
    input:
        os.path.join(out_dir, "Varscan", "VCF", "{sample_name}_varscan_snp.vcf"),
        os.path.join(out_dir, "Varscan", "VCF", "{sample_name}_varscan_indel.vcf")
    output:
        os.path.join(out_dir, "Varscan", "VCF", "{sample_name}_combined.vcf")
    params:
        python3 = config["python3"],
        out_dir = out_dir,
        base_path = base_path
    shell:
        """
        {params.python3} {params.base_path}/DNA/combine_vcf.py \
            --snv-vcf {params.out_dir}/Varscan/VCF/{wildcards.sample_name}_varscan_snp.vcf \
            --indel-vcf {params.out_dir}/Varscan/VCF/{wildcards.sample_name}_varscan_indel.vcf \
            --out-vcf {params.out_dir}/Varscan/VCF/{wildcards.sample_name}_combined.vcf
        """

# cmd13
rule check_step3:
    input:
        os.path.join(out_dir, "QC_Bam", "Bam", "{sample_name}_sort.mpileup"),
        os.path.join(out_dir, "Varscan", "VCF", "{sample_name}_varscan_snp.vcf"),
        os.path.join(out_dir, "Varscan", "VCF", "{sample_name}_varscan_indel.vcf"),
        os.path.join(out_dir, "Varscan", "VCF", "{sample_name}_combined.vcf")
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
            --workflow-name DNA \
            --step-num 3 \
            --process-log {params.out_dir}/TMP/{wildcards.sample_name}_Process.log && \
        touch {params.out_dir}/TMP/{wildcards.sample_name}_step3.done
        """