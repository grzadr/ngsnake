#!/bin/bash

-set eux

DOCKER_IMAGE = grzadr/biosak:mapping

docker pull ${DOCKER_IMAGE}
docker run -it \
  -v /etc/localtime:/etc/localtime:ro \
  --tmpfs /tmp:rw,exec,nosuid \
  -v ${1}:/data \
  -v ${PWD}:/home/jovyan/map \
  -v ${PWD}/Snakemake:/data/Snakemake \
  -w /home/jovyan/map \
  -it \
  --rm \
  ${DOCKER_IMAGE} /bin/bash
