#!/bin/bash

docker pull grzadr/biosak:mapping
docker run -it \
  -v /etc/localtime:/etc/localtime:ro \
  -v /data/OPUS:/data \
  -v ~/Git/Public/ngsake:/home/jovyan/map \
  -w /home/jovyan/map \
  -it \
  --rm \
  grzadr/biosak:mapping /bin/bash
