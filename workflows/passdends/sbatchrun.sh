#!/bin/bash

snakemake --cluster "sbatch -J snake" --jobs 15 --snakefile workflows/Snakefile_simdata --configfile workflows/passdends/config.yml --latency-wait 30
