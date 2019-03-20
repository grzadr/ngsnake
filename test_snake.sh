#!/bin/bash

set -eux

DOCKER_IMAGE=grzadr/biosak:OPUS2017SmallVariantCalling
NGSNAKE_DIR=${PWD}
SNPEFF_DIR=/data/SnpEff

#docker pull ${DOCKER_IMAGE}
docker run -it \
  -v /etc/localtime:/etc/localtime:ro \
  -v /tmp:/tmp:rw \
  -v ${1}:/data \
  -v ${NGSNAKE_DIR}:/ngsnake \
  -v ${SNPEFF_DIR}:/SnpEff \
  -w /data \
  -it \
  --rm \
  ${DOCKER_IMAGE} /bin/bash
