## Bionano - Clinical Affairs Simulation workflow
---
This pipeline is designed to simulate SVs based on a predefined file of SVs that adhere to the Simulation Format defined below. This workflow requires
Solve3.7.2 to be installed using conda/mamba in an environment named bionano_python3.0

## SETUP
---
Current implementation of this workflow utilizes the following packages
```
conda (4.12.0)
mamba (0.15.3)
Python (3.10.6)
Solve (3.7.2)
```

**Conda installation**

```
wget https://repo.anaconda.com/miniconda/Miniconda3-py39_4.12.0-Linux-x86_64.sh
#Follow-installation prompts
bash Miniconda3-py39_4.12.0-Linux-x86_64.sh
```

**Mamba installation**

```
conda install -n base -c conda-forge mamba
mamba create -c conda-forge -c bioconda -n snakemake snakemake
```
This will install snakemake into an isolated software environment, that has to be activated with

```
$ conda activate snakemake
$ snakemake --help
```

**Solve installation**

```
* Get Solve3.7.2
scp bionano@ansible:/home/bionano/Solve_build/Solve3.7_20221104_163244_27/solve_Solve3.7_r12432.12642rel_20221104_27.tar.gz . 

* Untar release file in installation directory
tar zxvf {solve.tar.gz}

* Create bionano_python3.0 conda/mamba environment
mamba env create -f envs/bionano_python3.0.yml

* Activate environment
conda activate bionano_python3.0 

* Change into installation directory
cd {installation directory}/bionano_packages (e.g. cd /usr/local/bioinformatics/common/releases/Solve3.7/Solve3.7_10192021_74_1/bionano_packages )

* Run the following commands
R CMD INSTALL ./BionanoR 
R CMD INSTALL ./FractCNV 
python3 -m pip install ./pybionano 
python3 -m pip install ./svconfmodels 
python3 -m pip install ./lohdetectiony 
```

**Install svconfmodels**
```
conda activate bionano_python3.0
git clone --branch develop git@bitbucket.org:biodiscovery/svconfmodels.git
cd svconfmodels
pip install .
```

**Install MolSim R Package from sv_simulation repository**
```
conda activate bionano_python3.0
git clone --branch develop git@bitbucket.org:biodiscovery/sv_simulation.git
cd sv_simulation/SingleMoleculeSimulator
R
install.packages("ghyp")
q()
R CMD INSTALL --no-multiarch --with-keep.source MolSim
```

***Installing into isolated environments is best practice in order to avoid compiler conflicts with non-dependent packages.***

----

## INTERPRETATION OF THE SIMULATION FORMAT: ##
Each SV is described by 7 values: (sv_type, acceptor chromosome, donor chromosome, breakpoint1, breakpoint2, size, orientation)
* chrom1 is always the acceptor chromosome/contig
* chrom2 is always the donor chromosome/contig
* All coordinates are EXCLUSIVE (the event happens after the start and before the end). This way insertions can be represented with the same convention as other events.
* Size is the actual size of the material deleted/translocated/inserted/duplicated/inverted


### deletions: ###
* breakpoint1: Start of the deletion (exclusive)
* breakpoint2 End of the deletion (exclusive)
* size: Ignored, inferred from the breakpoints: breakpoint2 - breakpoint1 - 1
* orientation: Ignored, always \+

### insertions: ###
* breakpoint1: Position immediately before the insertion (exclusive)
* breakpoint2: Position immediately after the insertion (exclusive): breakpoint1 + 1.
* size: The size of the insertion.
* ignored, always \+

### translocations: ###
* breakpoint1: Position in in the acceptor chromosome before the inserted material (exclusive)
* breakpoint2: Position in the donor chromosome before the the region translocated starts (exclusive)
* size: The size of the region to 'cut' from the donor
* Orientation: If +, the translocated region is kept as is. If -, the translocated region is inverted before 
  it is inserted in the acceptor position

### inversions: ###
* breakpoint1: Start of the inversion (exclusive)
* breakpoint2: End of the inversion (exclusive)
* size: Size of the inversion. It will be used when breakpoint2 is None
* Orientation: Ignored. The operation is done in the defined region as is

### duplications: ###
* breakpoint1: Position in in the acceptor chromosome that receives the duplicated region (exclusive,
    the inserted piece is added immediately after this point)
* breakpoint2: Position in in the donor chromosome where the duplicated region starts (exclusive)
* size: Size of the duplicated material.
* Orientation: If +, the region is attached as is. If -, the region attached is inverted first
* An example of tandem duplication of the region of 500bp chr1:1501-2000, added after position 2200: 
  (1,1, 2200, 1500, 500, \+)

----

**Example simulation format tsv**
----


| Type                   	| chrom1 	| chrom2 	| breakpoint1 	| breakpoint2 	| SVsize  	| Orientation 	|
|------------------------	|--------	|--------	|-------------	|-------------	|---------	|-------------	|
| translocation_intrachr 	| 1      	| 1      	| 3589895     	| 1E+08       	| 38580   	| +           	|
| duplication            	| 1      	| 1      	| 5250493     	| 5027224     	| 63322   	| +           	|
| deletion               	| 1      	| 1      	| 6682383     	| 6756792     	| 74408   	| +           	|
| inversion              	| 1      	| 1      	| 7337996     	| 7692969     	| 354972  	| +           	|
| inversion              	| 1      	| 1      	| 9342890     	| 9387508     	| 44617   	| +           	|
| translocation_interchr 	| 1      	| 10     	| 11237335    	| 11520342    	| 601140  	| +           	|
| duplication            	| 1      	| 1      	| 11431505    	| 11167103    	| 29346   	| +           	|
| translocation_interchr 	| 1      	| 15     	| 12487567    	| 94782974    	| 25164   	| +           	|
| translocation_intrachr 	| 1      	| 1      	| 14967516    	| 2.05E+08    	| 2617192 	| +           	|
| duplication            	| 1      	| 1      	| 15412855    	| 15133990    	| 231539  	| +           	|
| insertion              	| 1      	| 1      	| 15708791    	| 15708792    	| 617472  	| +           	|
| insertion              	| 1      	| 1      	| 15739110    	| 15739111    	| 1144353 	| +           	|
| insertion              	| 1      	| 1      	| 16950265    	| 16950266    	| 1082823 	| +           	|

--------
## Possible types of simulated structural variants

SV TYPE |
--------|
|cnv
|deletion
|deletion_nbase
|deletion_tiny
|duplication
|duplication_direct
|duplication_inverted
|duplication_split
|end
|gain
|insertion
|insertion_nbase
|insertion_tiny
|inversion
|inversion_nbase
|inversion_paired
|inversion_partial
|inversion_repeat
|loss
|monosomy
|translocation_interchr
|translocation_intrachr
|trans_interchr_common
|trans_intrachr_common
|trans_interchr_segdupe
|trans_intrachr_segdupe
|trans_interchr_repeat
|trans_intrachr_repeat
|trans_interchr_overlap
|trans_intrachr_overlap
|TP
|trisomy
|wild_type
|wild_type1
|wild_type2
|any
|complex
|rearrangements
|duplication_tandem

## RUNNING WORKFLOW
---

Do a dry-run of snakemake to ensure proper execution prior to launching the workflow.

```
screen
```

```
$(snakemake) snakemake -np --verbose
```
Execute snakemake workflow
```
snakemake -p --verbose --use-conda -c 1
```
Detach screen
```
ctrl-a ctrl-d
```