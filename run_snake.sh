#!/bin/bash

set -eux

DOCKER_IMAGE = grzadr/biosak:mapping

docker pull ${DOCKER_IMAGE}
docker run -it \
  -v /etc/localtime:/etc/localtime:ro \
  -v ${1}:/data \
  -v ${PWD}:/home/jovyan/map/ \
  -v ${PWD}/Snakemake:/data/Snakemake \
  -w /data \
  --name ngsnake_beta \
  --rm \
  --tmpfs /tmp:rw,exec,nosuid \
  ${DOCKER_IMAGE} \
  snakemake -s ~/map/Snakefile --configfile /data/config.yaml -pr -j ${2} "${@:3}"
