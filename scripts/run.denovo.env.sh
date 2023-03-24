#! /bin/bash
export SGE_ROOT=/home/sge
export DRMAA_LIBRARY_PATH=/home/sge/lib/lx-amd64/libdrmaa.so

if [[ $# -lt 2 ]]; then
	echo usage: $0 input.bnx output-dir [env]
	exit 1
fi

INPUT=$1
OUTPUT=$2
ENVIRONMENT=$3

ETC=/usr/local/bioinformatics/common/releases/Solve3.7/etc

if [[ -z $ENVIRONMENT ]]; then
	PIPELINE=/home/jestabrook/repositories/Solve3.7_20221104_163244_27
	CONDA=bionano_python3.0
	ANNOTATION=$ETC/annotation.hg38.ini
	REFERENCE=$PIPELINE/RefGenome/hg38_DLE1_0kb_0labels.cmap
else
	PIPELINE=/home/jestabrook/repositories/$ENVIRONMENT
	ANNOTATION=$ETC/annotation.${ENVIRONMENT}.ini
	REFERENCE=$PIPELINE/RefGenome/hg38_DLE1_0kb_0labels_masked_YPARs.cmap

	if [[ "$ENVIRONMENT" == "develop" ]]; then
		CONDA=bionano_dev_install
		REFERENCE=$PIPELINE/RefGenome/hg38_DLE1_0kb_0labels_masked_YPARs.cmap
	elif [[ "$ENVIRONMENT" == "Solve3.7" ]]; then
		CONDA=bionano_python3.0
		REFERENCE=$PIPELINE/RefGenome/hg38_DLE1_0kb_0labels.cmap
        elif [[ "$ENVIRONMENT" == "Solve3.7.1" ]]; then
		CONDA=Solve3.7.1
		REFERENCE=$PIPELINE/RefGenome/hg38_DLE1_0kb_0labels.cmap
	else
		echo unexpected environment: $ENVIRONMENT
		exit 1
	fi
fi

python3 $PIPELINE/Pipeline/1.0/pipelineCL.py \
  -l $OUTPUT \
  -t $PIPELINE/RefAligner/1.0 \
  -C $PIPELINE/Pipeline/1.0/clusterArguments_saphyr_supermicro.xml \
  -b $INPUT \
  -y \
  -d \
  -U \
  -i 5 \
  -F 1 \
  -W 1.0 \
  -c 0 \
  -a $PIPELINE/RefAligner/1.0/optArguments_haplotype_DLE1_saphyr_human.xml \
  -r $REFERENCE \
  -cm $PIPELINE/Process_Control_Datasets/FractCNV_Bed/data/homo_sapiens/hg38_cnv_masks.bed \
  -G $PIPELINE/Process_Control_Datasets/maskFiles/hg38_DLE1_gap_common_segdup_min10_com10kb_seg50kb.bed \
  --vapini $ANNOTATION \
  --compute-confidence human_hg38 \
  --species-reference human_hg38 \
  -f 0.08 \
  -J 64 \
  -Tn 48 \
  -j 64 \
  -jp 64 \
  -je 64 \
  -N 4 \
  --vap_conda $CONDA \
  --conda_source /home/jestabrook/miniconda3/etc/profile.d/conda.sh
