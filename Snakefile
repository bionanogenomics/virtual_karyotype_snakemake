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
        expand("samples/denovo/{sample}/pipeline_results.zip", sample=SAMPLES)

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
        # "samples/denovo/{sample}/data/annotation/variants_combine_filters_inMoleRefine1.smap",
        # "samples/denovo/{sample}/data/annotation/EXP_REFINEFINAL1.xmap",
        # "samples/denovo/{sample}/data/consensus_check/extension/sv/merged_smaps/EXP_REFINEFINAL1.cmap",
        # "samples/denovo/{sample}/data/auto_noise/autoNoise1_rescaled.bnx",
        # "samples/denovo/{sample}/data/auto_noise/autoNoise1.errbin",
        # "samples/denovo/{sample}/data/consensus_check/extension/sv/merged_smaps/EXP_REFINEFINAL1_q.cmap",
        # "samples/denovo/{sample}/data/consensus_check/extension/sv/merged_smaps/EXP_REFINEFINAL1_r.cmap",
        # "samples/denovo/{sample}/data/alignmolvref/copynumber/cnv_calls_exp.txt",
        # "samples/denovo/{sample}/data/alignmolvref/copynumber/cnv_chr_stats.txt"
    shell:
        """
        mkdir -p {params.out_dir}
        bash -c '
            . $HOME/.bashrc # if not loaded automatically
            conda activate bionano_python3.0
            ./scripts/run.denovo.env.sh {input} {params.out_dir} &> {log}
            conda deactivate'        
        """

