#!/bin/bash

docker pull grzadr/biosak:mapping
docker run -it \
  -v /etc/localtime:/etc/localtime:ro \
  -v ~/DSD:/data \
  -v ~/Git/Public/ngsnake/:/home/jovyan/map/ \
  -w /home/jovyan/map \
  -it \
  --rm \
  grzadr/biosak:mapping /bin/bash
