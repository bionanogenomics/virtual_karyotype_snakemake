__author__ = "Joey Estabrook"
__email__ = "jestabrook@bionano.com"

"""Bionano - Molecule simulation workflow"""

import os 
os.environ['OPENBLAS_NUM_THREADS'] = '1'
import numpy as np

configfile:"config.yaml"


def message(mes):
    sys.stderr.write("|--- " + mes + "\n")

dirs = ["scripts"]
for dir in dirs:
    if not os.path.exists(os.path.join(os.getcwd(),'logs',dir)):
        log_out = os.path.join(os.getcwd(), 'logs', dir)
        os.makedirs(log_out)

SAMPLES, = glob_wildcards('samples/{sample}/cnv_calls_exp.txt')

print(SAMPLES)

rule all:
    input:  
        expand("vk_results/{sample}/{sample}_results_ISCN.txt", sample=SAMPLES),
        expand("vk_results/{sample}/{sample}_custom_genome.txt", sample=SAMPLES),
        expand("vk_results/{sample}/{sample}_karyotype_rotated_split_4_of_4.pdf", sample=SAMPLES),
        expand("vk_results/{sample}/exp_refineFinal1_merged.xmap", sample=SAMPLES),
        "packaged_results/subset_postnatal_data.tar.gz",
        expand("vk_results/{sample}/{sample}_karyotype_rotated_merged_image.pdf", sample=SAMPLES)

rule run_vk:
    input:
        cnv="samples/{sample}/cnv_calls_exp.txt",
        smap="samples/{sample}/exp_refineFinal1_merged_filter_inversions.smap",
        rcmap="samples/{sample}/cnv_rcmap_exp.txt",
        xmap="samples/{sample}/exp_refineFinal1_merged.xmap"
    log:
        "logs/scripts/run_vk_{sample}.log"
    params:
        out_dir = lambda w: "vk_results/{sample}".format(sample=w.sample),
        centro = config['centro'],
        cyto = config['cytoband']
    conda:
        "envs/vk_env.yaml"
    resources:
        mem_mb=2000
    output:
        "vk_results/{sample}/{sample}_results.png",
        "vk_results/{sample}/{sample}_results.txt",
        "vk_results/{sample}/{sample}_results_ISCN.txt"
    shell:
        """
        python scripts/virtual_karyotype.py --cnv {input.cnv} --smap {input.smap} --rcmap {input.rcmap} --xmap {input.xmap} --centro {params.centro} --cyto {params.cyto} -n {wildcards.sample}_results -o {params.out_dir} &> {log}
        """

rule generate_plotting_parameter_files:
    input:
        "vk_results/{sample}/{sample}_results_ISCN.txt"
    log:
        "logs/scripts/generate_plotting_parameter_files_{sample}.log"
    conda:
        "envs/vk_env.yaml"
    params:
        cyto = config['resolved_cytoband']
    output:
        cytoband_out = "vk_results/{sample}/{sample}_custom_cytobands.txt",
        genome = "vk_results/{sample}/{sample}_custom_genome.txt",
        orientation = "vk_results/{sample}/{sample}_contig_orientation.txt",
        kprect = "vk_results/{sample}/{sample}_kprect_parameters.txt"
    shell:
        """
        python scripts/format_results.py --iscn_format {input} --cytoband {params.cyto} --cytobands_out {output.cytoband_out} --genome_out {output.genome} --contig_orientation_out {output.orientation} --kprect_out {output.kprect} &> {log}
        """

rule visualize_virtual_karyotype:
    input:
        cyto = "vk_results/{sample}/{sample}_custom_cytobands.txt",
        genome = "vk_results/{sample}/{sample}_custom_genome.txt",
        orientation = "vk_results/{sample}/{sample}_contig_orientation.txt",
        kprect = "vk_results/{sample}/{sample}_kprect_parameters.txt"
    log:
        "logs/scripts/visualize_virtual_karyotype{sample}.log"
    params:
        sample_handle = lambda w: "vk_results/{sample}/{sample}_".format(sample=w.sample),
    conda:
        "envs/karyoploter_env.yaml"
    output:
        "vk_results/{sample}/{sample}_processed.txt",
        temp("vk_results/{sample}/{sample}_karyotype_split_1_of_4.pdf"),
        temp("vk_results/{sample}/{sample}_karyotype_split_2_of_4.pdf"),
        temp("vk_results/{sample}/{sample}_karyotype_split_3_of_4.pdf"),
        temp("vk_results/{sample}/{sample}_karyotype_split_4_of_4.pdf")
    shell:
        """
        Rscript scripts/generate_ideogram.R {input.cyto} {input.genome} {input.orientation} {input.kprect} {output[0]} {params.sample_handle} &> {log}
        """

rule rotate_images:
    input:
        image1 = "vk_results/{sample}/{sample}_karyotype_split_1_of_4.pdf",
        image2 = "vk_results/{sample}/{sample}_karyotype_split_2_of_4.pdf",
        image3 = "vk_results/{sample}/{sample}_karyotype_split_3_of_4.pdf",
        image4 = "vk_results/{sample}/{sample}_karyotype_split_4_of_4.pdf"
    log:
        "logs/scripts/rotate_images_{sample}.log"
    conda:
        "envs/pypdf_env.yaml"
    output:
        "vk_results/{sample}/{sample}_karyotype_rotated_split_1_of_4.pdf",
        "vk_results/{sample}/{sample}_karyotype_rotated_split_2_of_4.pdf",
        "vk_results/{sample}/{sample}_karyotype_rotated_split_3_of_4.pdf",
        "vk_results/{sample}/{sample}_karyotype_rotated_split_4_of_4.pdf"
    shell:
        """
        python scripts/rotate_ideogram.py --image1 {input.image1} --image2 {input.image2} --image3 {input.image3} --image4 {input.image4} --out1 {output[0]} --out2 {output[1]} --out3 {output[2]} --out4 {output[3]} &> {log}
        """

rule merge_images:
    input:
        image1="vk_results/{sample}/{sample}_karyotype_rotated_split_1_of_4.pdf",
        image2="vk_results/{sample}/{sample}_karyotype_rotated_split_2_of_4.pdf",
        image3="vk_results/{sample}/{sample}_karyotype_rotated_split_3_of_4.pdf",
        image4="vk_results/{sample}/{sample}_karyotype_rotated_split_4_of_4.pdf"
    log:
        "logs/scripts/rotate_images_{sample}.log"
    conda:
        "envs/pypdf_env.yaml"
    output:
        "vk_results/{sample}/{sample}_karyotype_rotated_merged_image.pdf",
    shell:
        """
        python scripts/merge_images.py --image1 {input.image1} --image2 {input.image2} --image3 {input.image3} --image4 {input.image4} --out1 {output[0]} &> {log}
        """

rule copy_data:
    input:
        cnv="samples/{sample}/cnv_calls_exp.txt",
        smap="samples/{sample}/exp_refineFinal1_merged_filter_inversions.smap",
        rcmap="samples/{sample}/cnv_rcmap_exp.txt",
        xmap="samples/{sample}/exp_refineFinal1_merged.xmap"
    log:
        "logs/scripts/copy_data_{sample}.log"
    params:
        out_dir = lambda w: "vk_results/{sample}".format(sample=w.sample),
    output:
        "vk_results/{sample}/cnv_calls_exp.txt",
        "vk_results/{sample}/exp_refineFinal1_merged_filter_inversions.smap",
        "vk_results/{sample}/cnv_rcmap_exp.txt",
        "vk_results/{sample}/exp_refineFinal1_merged.xmap"
    shell:
        """
        cp {input.cnv} {params.out_dir}
        cp {input.smap} {params.out_dir}
        cp {input.rcmap} {params.out_dir}
        cp {input.xmap} {params.out_dir}
        """

rule package_results:
    input:
        "vk_results",
    log:
        "logs/scripts/package_results.log"
    output:
        data_out="packaged_results/subset_postnatal_data.tar.gz",
    shell:
        """
        mkdir -p packaged_results
        tar -zcvf {output.data_out} {input}
        """