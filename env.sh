#!/bin/bash

# Export WESTPA environment on Cluster
export WEST_ROOT=/projects/lemkul_lab/software/tinkercliffs/openmm/7.7.0/conda_envs/openmm7.7.0

# Activate the openmm environment on computing nodes
module reset 
module load Anaconda3/2022.05
echo "Checking Anaconda loading..."
which python
source activate /projects/lemkul_lab/software/tinkercliffs/openmm/7.7.0/conda_envs/openmm7.7.0
echo "Environment should be loaded..."
which python

export WM_ZMQ_MASTER_HEARTBEAT=100
export WM_ZMQ_WORKER_HEARTBEAT=100
export WM_ZMQ_TIMEOUT_FACTOR=200
export WEST_SIM_ROOT="$PWD"
export SIM_NAME=$(basename $WEST_SIM_ROOT)
export BASH=/bin/bash
export PERL=/usr/bin/perl
export ZSH=/bin/zsh
export IFCONFIG=/bin/ifconfig
export CUT=/usr/bin/cut
export TR=/usr/bin/tr
export LN=/bin/ln
export CP=/bin/cp
export RM=/bin/rm
export SED=/bin/sed
export CAT=/bin/cat
export HEAD=/bin/head
export TAR=/bin/tar
export AWK=/usr/bin/awk
export PASTE=/usr/bin/paste
export GREP=/bin/grep
export SORT=/usr/bin/sort
export UNIQ=/usr/bin/uniq
export HEAD=/usr/bin/head
export MKDIR=/bin/mkdir
export ECHO=/bin/echo
export DATE=/bin/date
