import os
out_dir = config["out_dir"]

rule fusion_anno:
    input:
        reads_stat = os.path.join(out_dir, "TMP", "{sample_name}_reads.stats")
    output:
        os.path.join(out_dir, "TMP", "{sample_name}_fusion_raw.txt")
    params:
        python3 = config["python3"],
        base_path = config["base_path"],
        fusion_anno_txt = config['RNA']["fusion_anno_txt"]
    shell:
        """
        {params.python3} {params.base_path}/RNA/fusion_anno.py \
			--fusion-reads {input.reads_stat} \
			--fusion-anno {params.fusion_anno_txt} \
			--sample-name {wildcards.sample_name} \
			--out-file {output}
        """


