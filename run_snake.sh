#!/bin/bash

docker pull grzadr/biosak:mapping
docker run -it \
  -v /etc/localtime:/etc/localtime:ro \
  -v ${1}:/data \
  -v ${PWD}:/home/jovyan/map \
  -w /home/jovyan/map \
  --name ngsnake \
  grzadr/biosak:mapping snakemake -pr -j ${2} --configfile ${3} ${4}
