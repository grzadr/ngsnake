#!/bin/bash

set -eux

docker run -it \
  --pull \
  -v /etc/localtime:/etc/localtime:ro \
  -v ${1}:/data \
  -v ${PWD}:/home/jovyan/map/ \
  -w /data \
  --name ngsnake_beta \
  --rm \
  --tmpfs /tmp:rw,exec,nosuid \
  grzadr/biosak:mapping \
  snakemake -s ~/map/Snakefile --configfile /data/config.yaml -pr -j ${2} "${@:3}"
