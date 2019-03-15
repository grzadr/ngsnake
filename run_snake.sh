#!/bin/bash

set -eux

DOCKER_IMAGE=grzadr/biosak:OPUS2017SmallVariantCalling

THREADS=20
DATA_DIR=${1}
SNAKEFILE=${2}
NGSNAKE_DIR=${PWD}
SNPEFF_DIR=/data/SnpEff
SNAKEMAKE_ARGS="${@:3}"

#docker pull ${DOCKER_IMAGE}
docker run -it \
  -v /etc/localtime:/etc/localtime:ro \
  -v ${DATA_DIR}:/data \
  -v ${NGSNAKE_DIR}:/ngsnake \
  -v ${SNPEFF_DIR}:/SnpEff \
  -w /data \
  --name ngsnake_mapping \
  --rm \
  -v /tmp:/tmp:rw \
  ${DOCKER_IMAGE} \
  snakemake \
  -s "/ngsnake/${SNAKEFILE}" \
  --configfile /ngsnake/config.yaml \
  --resources mem_mb=188416 \
  -pr -j ${THREADS} \
  ${SNAKEMAKE_ARGS}

