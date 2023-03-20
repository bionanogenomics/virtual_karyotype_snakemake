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
            touch {output[0]}
            mv {params.out_dir}/*.bnx {output[1]}
            conda deactivate'
        """