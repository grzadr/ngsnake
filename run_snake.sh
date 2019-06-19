#!/bin/bash

set -eux

DOCKER_IMAGE_TAG="${1}"
DATA_DIR=${2}
SNAKE_FILE=${3}
SNAKE_CONFIG=${4}
SNAKE_ARGS="${@:5}"

DOCKER_IMAGE="grzadr/biosak:${DOCKER_IMAGE_TAG}"
DOCKER_CONTAINER="ngsnake_$(date '+%Y-%m-%d_%H%M%S')"
SNAKE_THREADS=20
SNAKE_MEMORY=188416
SNAKE_DIR=${PWD}
SNPEFF_DIR=/data/SnpEff

#docker pull ${DOCKER_IMAGE}
docker run -it \
  -v /etc/localtime:/etc/localtime:ro \
  -v ${DATA_DIR}:/data \
  -v ${SNAKE_DIR}:/ngsnake \
  -v ${SNPEFF_DIR}:/SnpEff \
  -w /data \
  --name ${DOCKER_CONTAINER} \
  --rm \
  -v /tmp:/tmp:rw \
  ${DOCKER_IMAGE} \
  snakemake \
  -s "${SNAKE_FILE}" \
  --configfile "${SNAKE_CONFIG}" \
  --resources mem_mb=${SNAKE_MEMORY} \
  -pr -j ${SNAKE_THREADS} \
  ${SNAKE_ARGS}

