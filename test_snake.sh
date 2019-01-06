#!/bin/bash

set -eux

DOCKER_IMAGE=grzadr/biosak:OPUS-2017-DSD

#docker pull ${DOCKER_IMAGE}
docker run -it \
  -v /etc/localtime:/etc/localtime:ro \
  --tmpfs /tmp:rw,exec,nosuid \
  -v ${1}:/data \
  -v ${PWD}:/ngsnake \
  -w /data \
  -it \
  --rm \
  ${DOCKER_IMAGE} /bin/bash
