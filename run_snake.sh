#!/bin/bash

set -eux

DOCKER_IMAGE=grzadr/biosak:OPUS-2017-DSD_variant_calling

THREADS=20
DATA_DIR=${1}
SNAKEMAKE_ARGS="${@:2}"

docker pull ${DOCKER_IMAGE}
docker run -it \
  -v /etc/localtime:/etc/localtime:ro \
  -v ${DATA_DIR}:/data \
  -v ${PWD}:/ngsnake \
  -w /data \
  --name ngsnake_mapping \
  --rm \
  --tmpfs /tmp:rw,exec,nosuid \
  ${DOCKER_IMAGE} \
  snakemake -s /ngsnake/Snakefile --resources mem_mb=188416 --configfile /ngsnake/config.yaml -pr -j ${THREADS} ${SNAKEMAKE_ARGS}
