import os

rule clear_result:
    input:
        os.path.join(out_dir, "QC_Bam", "Bam", "{sample_name}.bam"),
        os.path.join(out_dir, "Run_scripts", "{sample_name}.sh")
    output:
        os.path.join(out_dir, "QC_Bam", "Bam", "clear.done")
    params:
        python3 = config["python3"],
        base_path = config["base_path"],
        out_dir = config["out_dir"]
    shell:
        """
        {params.python3}  {params.base_path}/common/clear_file.py --file-list {input} && \
            touch {params.out_dir}/QC_Bam/Bam/clear.done
        """


