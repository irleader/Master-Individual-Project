#!/bin/bash
# Copyright 2018 Google LLC.
#modified by Jiajia Xu.2019.

#set or unset options and positional parameters,
#-e When this option is on, if a simple command fails for any of the reasons listed in Consequences of Shell  Errors or returns an exit status value >0, and is not part of the compound list following a while, until, or if keyword, and is not a part of an AND  or  OR  list,  and  is  not  a  pipeline preceded by the ! reserved word, then the shell shall immediately exit.
#-u The shell shall write a message to standard error when it tries to expand a variable that  is  not set and immediately exit. An interactive shell shall not exit.
#-o     Write the current settings of the options to standard output in an unspecified format.

set -euo pipefail

## Preliminaries
# Set a number of shell variables, to make what follows easier to read.
#base directory
BASE="${PWD}/wgs-case-study"
#deepvariant version no.
BIN_VERSION="0.8.0"
#input directory
INPUT_DIR="${BASE}/input/data"
#reference genome in fasta gunzip format
REF="hs37d5.fa.gz"
#aligned reads in bam format
BAM="HG002_NIST_150bp_50x.bam"
#True variants
TRUTH_VCF="HG002_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-22_v.3.3.2_highconf_triophased.vcf.gz"
TRUTH_BED="HG002_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-22_v.3.3.2_highconf_noinconsistent.bed"
#No of CPUs
N_SHARDS="12"
#output directory
OUTPUT_DIR="${BASE}/output"
#called variants
OUTPUT_VCF="HG002.output.vcf.gz"
OUTPUT_GVCF="HG002.output.g.vcf.gz"
LOG_DIR="${OUTPUT_DIR}/logs"

## Create local directory structure
mkdir -p "${OUTPUT_DIR}"
mkdir -p "${INPUT_DIR}"
mkdir -p "${LOG_DIR}"


## Download extra packages
# Install aria2 to download data files.
sudo apt-get -qq -y update
sudo apt-get -qq -y install aria2

# if already have docker, can comment out
if ! hash docker 2>/dev/null; then
  echo "'docker' was not found in PATH. Installing docker..."
  # Install docker using instructions on:
  # https://docs.docker.com/install/linux/docker-ce/ubuntu/
  sudo apt-get -qq -y install \
    apt-transport-https \
    ca-certificates \
    curl \
    gnupg-agent \
    software-properties-common
  curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo apt-key add -
  sudo add-apt-repository \
    "deb [arch=amd64] https://download.docker.com/linux/ubuntu \
    $(lsb_release -cs) \
    stable"
  sudo apt-get -qq -y update
  sudo apt-get -qq -y install docker-ce
fi


# Copy the data
aria2c -c -x10 -s10 "https://storage.googleapis.com/deepvariant/case-study-testdata/${TRUTH_BED}" -d "${INPUT_DIR}"
aria2c -c -x10 -s10 "https://storage.googleapis.com/deepvariant/case-study-testdata/${TRUTH_VCF}" -d "${INPUT_DIR}"
aria2c -c -x10 -s10 "https://storage.googleapis.com/deepvariant/case-study-testdata/${TRUTH_VCF}.tbi" -d "${INPUT_DIR}"
aria2c -c -x10 -s10 "https://storage.googleapis.com/deepvariant/case-study-testdata/${BAM}" -d "${INPUT_DIR}"
aria2c -c -x10 -s10 "https://storage.googleapis.com/deepvariant/case-study-testdata/${BAM}.bai" -d "${INPUT_DIR}"
aria2c -c -x10 -s10 https://storage.googleapis.com/deepvariant/case-study-testdata/hs37d5.fa.gz -d "${INPUT_DIR}"
aria2c -c -x10 -s10 https://storage.googleapis.com/deepvariant/case-study-testdata/hs37d5.fa.gz.fai -d "${INPUT_DIR}"
aria2c -c -x10 -s10 https://storage.googleapis.com/deepvariant/case-study-testdata/hs37d5.fa.gz.gzi -d "${INPUT_DIR}"
aria2c -c -x10 -s10 https://storage.googleapis.com/deepvariant/case-study-testdata/hs37d5.fa.gzi -d "${INPUT_DIR}"
aria2c -c -x10 -s10 https://storage.googleapis.com/deepvariant/case-study-testdata/hs37d5.fa.fai -d "${INPUT_DIR}"

#already have deepvariant, can comment out
# Pull the docker image.
sudo docker pull gcr.io/deepvariant-docker/deepvariant:"${BIN_VERSION}"

echo "Run DeepVariant..."
sudo docker run \
  -v "${INPUT_DIR}":"/input" \
  -v "${OUTPUT_DIR}:/output" \
  gcr.io/deepvariant-docker/deepvariant:"${BIN_VERSION}" \
  /opt/deepvariant/bin/run_deepvariant \
  --model_type=WGS \
  --ref="/input/${REF}" \
  --reads="/input/${BAM}" \
  --output_vcf=/output/${OUTPUT_VCF} \
  --output_gvcf=/output/${OUTPUT_GVCF} \
  --num_shards=${N_SHARDS}
echo "Done."
echo

## Evaluation: run hap.py
echo "Start evaluation with hap.py..."
UNCOMPRESSED_REF="${INPUT_DIR}/hs37d5.fa"

# hap.py cannot read the compressed fa, so uncompress
# into a writable directory. Index file was downloaded earlier.
zcat <"${INPUT_DIR}/${REF}" >"${UNCOMPRESSED_REF}"

sudo docker pull pkrusche/hap.py
( sudo docker run -i \
-v "${INPUT_DIR}:${INPUT_DIR}" \
-v "${OUTPUT_DIR}:${OUTPUT_DIR}" \
pkrusche/hap.py /opt/hap.py/bin/hap.py \
  "${INPUT_DIR}/${TRUTH_VCF}" \
  "${OUTPUT_DIR}/${OUTPUT_VCF}" \
  -f "${INPUT_DIR}/${TRUTH_BED}" \
  -r "${UNCOMPRESSED_REF}" \
  -l 20 \
  -o "${OUTPUT_DIR}/happy20.output" \
  --engine=vcfeval
) 2>&1 | tee "${LOG_DIR}/happy20.log"
echo "Done."
