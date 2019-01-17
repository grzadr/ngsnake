#!/bin/bash

set -eux

DOCKER_IMAGE=grzadr/biosak:mapping

docker pull ${DOCKER_IMAGE}
docker run -it \
  -v /etc/localtime:/etc/localtime:ro \
  -v /tmp:/tmp:rw \
  -v ${1}:/data \
  -v ${PWD}:/ngsnake \
  -w /data \
  -it \
  --rm \
  ${DOCKER_IMAGE} /bin/bash
