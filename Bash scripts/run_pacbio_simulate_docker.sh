#!/bin/bash
# Copyright 2018 Google LLC.
#modified by Jiajia Xu.2019.

#set or unset options and positional parameters,
#-e When this option is on, if a simple command fails for any of the reasons listed in Consequences of Shell  Errors or returns an exit status value >0, and is not part of the compound list following a while, until, or if keyword, and is not a part of an AND  or  OR  list,  and  is  not  a  pipeline preceded by the ! reserved word, then the shell shall immediately exit.
#-u The shell shall write a message to standard error when it tries to expand a variable that  is  not set and immediately exit. An interactive shell shall not exit.
#-o     Write the current settings of the options to standard output in an unspecified format.

set -eu pipefail

## Preliminaries
# Set a number of shell variables, to make what follows easier to read.

BASE="${PWD}/pacbio_simulate"
BIN_VERSION="0.8.0"
INPUT_DIR="${BASE}/input/data"
REF="chr20_hs37d5.fa.gz"
BAM="sorted.p.SNP_INDEL_chr201.bam"
TRUTH_VCF="sorted_SNP_INDEL_chr1_22_hs37d5.vcf.gz"
#TRUTH_BED="HG002_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-22_v.3.3.2_highconf_noinconsistent.bed"

N_SHARDS="12"
OUTPUT_DIR="${BASE}/output"
OUTPUT_VCF="p.SNP_INDEL_chr201.output.vcf.gz"
OUTPUT_GVCF="p.SNP_INDEL_chr201.output.g.vcf.gz"
LOG_DIR="${OUTPUT_DIR}/logs"

mkdir -p "${OUTPUT_DIR}"
mkdir -p "${INPUT_DIR}"
mkdir -p "${LOG_DIR}"


# Pull the docker image.
sudo docker pull gcr.io/deepvariant-docker/deepvariant:"${BIN_VERSION}"

echo "Run DeepVariant..."
sudo docker run \
  -v "${INPUT_DIR}":"/input" \
  -v "${OUTPUT_DIR}:/output" \
  gcr.io/deepvariant-docker/deepvariant:"${BIN_VERSION}" \
  /opt/deepvariant/bin/run_deepvariant \
  --model_type=PACBIO \
  --ref="/input/${REF}" \
  --reads="/input/${BAM}" \
  --output_vcf=/output/${OUTPUT_VCF} \
  --output_gvcf=/output/${OUTPUT_GVCF} \
  --num_shards=${N_SHARDS}
echo "Done."
echo

## Evaluation: run hap.py
echo "Start evaluation with hap.py..."
UNCOMPRESSED_REF="${INPUT_DIR}/chr20_hs37d5.fa"

# hap.py cannot read the compressed fa, so uncompress
# into a writable directory. Index file was downloaded earlier.
#zcat <"${INPUT_DIR}/${REF}" >"${UNCOMPRESSED_REF}"

sudo docker pull pkrusche/hap.py
( sudo docker run -i \
-v "${INPUT_DIR}:${INPUT_DIR}" \
-v "${OUTPUT_DIR}:${OUTPUT_DIR}" \
pkrusche/hap.py /opt/hap.py/bin/hap.py \
"${INPUT_DIR}/${TRUTH_VCF}" \
"${OUTPUT_DIR}/${OUTPUT_VCF}" \
-r "${UNCOMPRESSED_REF}" \
-l 20 \
-o "${OUTPUT_DIR}/p.SNP_INDEL_chr201.happy.output" \
--engine=vcfeval
) 2>&1 | tee "${LOG_DIR}/p.SNP_INDEL_chr201.happy.log"
echo "Done."
