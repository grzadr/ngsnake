#!/bin/bash

docker pull grzadr/biosak:mapping
docker run -it \
  -v /etc/localtime:/etc/localtime:ro \
  -v ${1}:/data \
  -v ${PWD}/Snakefile:/home/jovyan/map/ \
  -w /home/jovyan/map \
  --name ngsnake \
  grzadr/biosak:mapping snakemake --configfile /data/config.yaml -pr -j ${2} ${3}
