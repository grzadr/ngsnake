#!/bin/bash

set -eux

DOCKER_IMAGE = grzadr/biosak:mapping

docker pull ${DOCKER_IMAGE}
docker run -it \
  -v /etc/localtime:/etc/localtime:ro \
  -v ${1}:/data \
  -v ${PWD}:/ngsnake \
  -w /data \
  --name ngsnake_mapping \
  --rm \
  --tmpfs /tmp:rw,exec,nosuid \
  ${DOCKER_IMAGE} \
  snakemake -s /ngsnake/Snakefile --resources mem_mb=196608 --configfile /data/config.yaml -pr -j ${2} "${@:3}"
