#!/bin/bash

cp ${PWD}/Snakefile ${1}/Snakefile
docker pull grzadr/biosak:mapping
docker run -it \
  -v /etc/localtime:/etc/localtime:ro \
  -v ${1}:/data \
  -v ${PWD}/Snakefile:/home/jovyan/Snakefile \
  -it \
  --rm \
  grzadr/biosak:mapping /bin/bash
