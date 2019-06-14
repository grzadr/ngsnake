#!/bin/bash

DATADIR="/data/OPUS-2017-Obesity"

./run_snake.sh OPUS-Obesity-2017-SmallVariantCalling "$DATADIR" /ngsnake/Snakefile.small_variants /ngsnake/config.yaml multiqc_mapping_report
./run_snake.sh OPUS-Obesity-2017-SmallVariantCalling "$DATADIR" /ngsnake/Snakefile.small_variants /ngsnake/config.yaml call_variants
./run_snake.sh OPUS-Obesity-2017-SmallVariantCalling "$DATADIR" /ngsnake/Snakefile.copy_number_variants /ngsnake/config.yaml call_cnvs
./run_snake.sh OPUS-Obesity-2017-SmallVariantCalling "$DATADIR" /ngsnake/Snakefile.copy_number_variants /ngsnake/config.yaml call_cnmops_cvsc

