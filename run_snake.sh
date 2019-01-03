#!/bin/bash

set -eux

# echo "${@:3}"
# exit 0
docker pull grzadr/biosak:mapping
docker run -it \
  -v /etc/localtime:/etc/localtime:ro \
  -v ${1}:/data \
  -v ${PWD}:/home/jovyan/map/ \
  -w /home/jovyan/map \
  --name ngsnake_beta \
  --rm \
  --tmpfs /tmp:rw,exec,nosuid \
  grzadr/biosak:mapping snakemake --configfile /data/config.yaml -pr -j ${2} "${@:3}"
