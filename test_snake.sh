#!/bin/bash

-set eux

docker run -it \
  --pull \
  -v /etc/localtime:/etc/localtime:ro \
  --tmpfs /tmp:rw,exec,nosuid \
  -v ${1}:/data \
  -v ${PWD}:/home/jovyan/map \
  -w /home/jovyan/map \
  -it \
  --rm \
  grzadr/biosak:mapping /bin/bash
