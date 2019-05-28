#!/bin/bash

set -eux

DOCKER_IMAGE_TAG="${1}"
DATA_DIR=${2}

DOCKER_IMAGE="grzadr/biosak:${DOCKER_IMAGE_TAG}"
SNAKE_DIR=${PWD}
SNPEFF_DIR=/data/SnpEff

#docker pull ${DOCKER_IMAGE}
docker run -it \
  -v /etc/localtime:/etc/localtime:ro \
  -v /tmp:/tmp:rw \
  -v ${DATA_DIR}:/data \
  -v ${SNAKE_DIR}:/ngsnake \
  -v ${SNPEFF_DIR}:/SnpEff \
  -w /data \
  -it \
  --rm \
  ${DOCKER_IMAGE} /bin/bash
