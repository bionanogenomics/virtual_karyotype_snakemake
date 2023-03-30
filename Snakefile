__author__ = "Joey Estabrook"
__email__ = "jestabrook@bionano.com"

"""Bionano - Molecule simulation workflow"""

import os 

configfile:"config.yaml"

dirs = ["scripts"]
print("Generating log directories....")
for dir in dirs:
    if not os.path.exists(os.path.join(os.getcwd(),'logs',dir)):
        log_out = os.path.join(os.getcwd(), 'logs', dir)
        os.makedirs(log_out)
        print(log_out)
print("Finished!")

SAMPLES = [1]

rule all:
    input:  
        expand("samples/{sample}/molecules/task_complete.tsv", sample=SAMPLES),
        expand("samples/denovo/{sample}/pipeline_results.zip", sample=SAMPLES),
        expand("vk_results/{sample}/{sample}_results_ISCN.txt", sample=SAMPLES),
        expand("vk_results/{sample}/{sample}_custom_genome.txt", sample=SAMPLES),
        expand("vk_results/{sample}/{sample}_karyotype_rotated_split_1_of_2.pdf", sample=SAMPLES),
        expand("vk_results/{sample}/{sample}_karyotype_rotated_split_2_of_2.pdf", sample=SAMPLES)

rule select_svs:
    input:
        "input/hg19_random_svs_2.txt",
    log:
        "logs/scripts/select_svs_{sample}.log"
    params:
        out_dir = lambda w: "samples/{sample}".format(sample=w.sample),
    output:
        "samples/{sample}/{sample}_sv_truth_set.tsv"
    shell:
        """
        mkdir -p {params.out_dir}
        python scripts/sv_selector.py --input {input[0]} --output {output[0]} &> {log}
        """

rule generate_cmap:
    input:
        "samples/{sample}/{sample}_sv_truth_set.tsv"
    log:
        "logs/scripts/generat_cmap_{sample}.log"
    params:
        reference_cmap = config['reference_cmap'],
        bngomodel = config['bngomodel']
    output:
        "samples/{sample}/{sample}_simulation.cmap",
    shell:
        """
        bash -c '
            . $HOME/.bashrc # if not loaded automatically
            conda activate bionano_python3.0
            {params.bngomodel} svs apply --genome {params.reference_cmap} --svs {input[0]} --output {output[0]} &> {log}
            conda deactivate' """

rule simulate_molecules:
    input:
        cmap = "samples/{sample}/{sample}_simulation.cmap",
    log:
        "logs/scripts/simulate_molecules_{sample}.log"
    params:
        out_dir = lambda w: "samples/{sample}/molecules".format(sample=w.sample),
        sim_params = config['sim_params']
    output:
        "samples/{sample}/molecules/task_complete.tsv",
        "samples/{sample}/molecules/{sample}_molecules.bnx"
    shell:
        """
        mkdir -p {params.out_dir}
        bash -c '
            . $HOME/.bashrc # if not loaded automatically
            conda activate bionano_python3.0
            Rscript scripts/simulate_molecules.R -r {input.cmap} -o {params.out_dir} -p {params.sim_params} &> {log}
            mv {params.out_dir}/*mresSD*.bnx {output[1]}
            touch {output[0]}
            conda deactivate'
        """

rule run_denovo:
    input:
        "samples/{sample}/molecules/{sample}_molecules.bnx"
    log:
        "logs/scripts/run_denovo_{sample}.log"
    params:
        out_dir = lambda w: "samples/denovo/{sample}".format(sample=w.sample),
    output:
        "samples/denovo/{sample}/pipeline_results.zip",
        "samples/denovo/{sample}/data/annotation/variants_combine_filters_inMoleRefine1.smap",
        "samples/denovo/{sample}/data/annotation/EXP_REFINEFINAL1.xmap",
        "samples/denovo/{sample}/data/consensus_check/extension/sv/merged_smaps/EXP_REFINEFINAL1.cmap",
        "samples/denovo/{sample}/data/auto_noise/autoNoise1_rescaled.bnx",
        "samples/denovo/{sample}/data/auto_noise/autoNoise1.errbin",
        "samples/denovo/{sample}/data/consensus_check/extension/sv/merged_smaps/EXP_REFINEFINAL1_q.cmap",
        "samples/denovo/{sample}/data/consensus_check/extension/sv/merged_smaps/EXP_REFINEFINAL1_r.cmap",
        "samples/denovo/{sample}/data/alignmolvref/copynumber/cnv_calls_exp.txt",
        "samples/denovo/{sample}/data/alignmolvref/copynumber/cnv_chr_stats.txt"
    shell:
        """
        mkdir -p {params.out_dir}
        bash -c '
            . $HOME/.bashrc # if not loaded automatically
            conda activate bionano_python3.0
            ./scripts/run.denovo.env.sh {input} {params.out_dir} &> {log}
            conda deactivate'        
        """

rule run_vk:
    input:
        cnv="samples/denovo/{sample}/contigs/alignmolvref/copynumber/cnv_calls_exp.txt",
        smap="samples/denovo/{sample}/contigs/exp_refineFinal1_sv/merged_smaps/exp_refineFinal1_merged_filter_inversions.smap",
        rcmap="samples/denovo/{sample}/contigs/alignmolvref/copynumber/cnv_rcmap_exp.txt",
        xmap="samples/denovo/{sample}/contigs/annotation/exp_refineFinal1_merged.xmap"
    log:
        "logs/scripts/run_vk_{sample}.log"
    params:
        out_dir = lambda w: "vk_results/{sample}".format(sample=w.sample),
        centro = config['centro'],
        cyto = config['cytoband']
    conda:
        "envs/vk_env.yaml"
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
        cytoband_out = "vk_results/{sample}/{sample}_custom_cytobands.txt",
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
        temp("vk_results/{sample}/{sample}_karyotype_split_1_of_2.pdf"),
        temp("vk_results/{sample}/{sample}_karyotype_split_2_of_2.pdf")
    shell:
        """
        Rscript scripts/generate_ideogram.R {input.cytoband_out} {input.genome} {input.orientation} {input.kprect} {output[0]} {params.sample_handle} &> {log}
        """

rule rotate_images:
    input:
        image1 = "vk_results/{sample}/{sample}_karyotype_split_1_of_2.pdf",
        image2 = "vk_results/{sample}/{sample}_karyotype_split_2_of_2.pdf"
    log:
        "logs/scripts/rotate_images_{sample}.log"
    conda:
        "envs/pypdf_env.yaml"
    output:
        "vk_results/{sample}/{sample}_karyotype_rotated_split_1_of_2.pdf",
        "vk_results/{sample}/{sample}_karyotype_rotated_split_2_of_2.pdf"
    shell:
        """
        python scripts/rotate_ideogram.py --image1 {input.image1} --image2 {input.image2} --out1 {output[0]} --out2 {output[1]} &> {log}
        """