#!/bin/bash

# Make sure environment is set
source env.sh

# Clean up
rm -f west.log
#conda activate westpa-2.0
#source activate /projects/lemkul_lab/software/tinkercliffs/openmm/7.7.0/conda_envs/openmm7.7.0

# Run w_run
#w_run --work-manager processes "$@" &> west.log
w_run --work-manager serial > west.log
