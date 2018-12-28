#!/bin/bash

#cp ${PWD}/Snakefile ${1}/Snakefile
docker pull grzadr/biosak:mapping
docker run -it \
  -v /etc/localtime:/etc/localtime:ro \
  -v ${1}:/data \
  -v ${PWD}:/home/jovyan/map \
  -w /home/jovyan/map \
  -it \
  --rm \
  grzadr/biosak:mapping /bin/bash
