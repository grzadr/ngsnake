#!/bin/bash

set -eux

DOCKER_IMAGE=grzadr/biosak:mapping

THREADS=20
DATA_DIR=${1}
NGSNAKE_DIR=${PWD}
SNPEFF_DIR=/data/SnpEff
SNAKEMAKE_ARGS="${@:2}"

docker pull ${DOCKER_IMAGE}
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
  snakemake -s /ngsnake/Snakefile --resources mem_mb=188416 --configfile /ngsnake/config.yaml -pr -j ${THREADS} ${SNAKEMAKE_ARGS}

