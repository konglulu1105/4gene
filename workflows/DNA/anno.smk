import os
out_dir = config["out_dir"]
base_path = config["base_path"]
python3 = config["python3"]

rule anno_snv:
    input:
        os.path.join(out_dir, "Varscan", "VCF", "{sample_name}_combined.vcf")
    output:
        os.path.join(out_dir, "TMP", "{sample_name}_SNV_anno.xls")
    params:
        python3 = python3,
        base_path = base_path,
        out_dir = out_dir,
        zc_hotspots = config["DNA"]["zc_hotspots"]
    shell:
        """
        {params.python3} {params.base_path}/DNA/anno_snv.py \
            --hotspot {params.zc_hotspots} \
            --vcf {input} \
            --sample-name {wildcards.sample_name} \
            --out-file {params.out_dir}/TMP/{wildcards.sample_name}_SNV_anno.xls
        """

rule add_extra_snv:
    input:
        os.path.join(out_dir, "TMP", "{sample_name}_SNV_anno.xls")
    output:
        os.path.join(out_dir, "results", "{sample_name}_SNV.xls")
    params:
        python3 = python3,
        base_path = base_path,
        out_dir = out_dir
    shell:
        """
        {params.python3} {params.base_path}/common/add_extra_sites.py \
            --inp-file {params.out_dir}/TMP/{wildcards.sample_name}_SNV_anno.xls \
            --out-file {params.out_dir}/results/{wildcards.sample_name}_SNV.xls \
            --variant-type SNV \
            --sample-name {wildcards.sample_name}
        """
