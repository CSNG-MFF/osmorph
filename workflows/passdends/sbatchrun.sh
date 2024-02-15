#!/bin/bash

snakemake --cluster "sbatch -J snake -w w12" --jobs 16 --snakefile workflows/Snakefile_simdata --configfile workflows/passdends/config.yml --latency-wait 30
