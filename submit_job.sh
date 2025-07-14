#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --gres=gpu:8
#SBATCH -t 144:00:00
#SBATCH -p dgx_normal_q
#SBATCH -A lemkulgroup
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mdpoleto@vt.edu
#SBATCH --job-name=abl1
#SBATCH --cpus-per-task=4
#SBATCH --exclusive


module reset # <<<< MUST BE RESET, NOT PURGE
module load Anaconda3/2022.05

echo ">>> Anaconda should be loaded:"
which python

# Activate my Conda Enviroment with openMM+westpa
source activate /projects/lemkul_lab/software/tinkercliffs/openmm/7.7.0/conda_envs/openmm7.7.0

echo ">>> Environment should be loaded:"
which python

echo "################################################"
##################################
set -x
cd $SLURM_SUBMIT_DIR
source env.sh || exit 1

env | sort

cd $WEST_SIM_ROOT


SERVER_INFO=$WEST_SIM_ROOT/west_zmq_info-$SLURM_JOBID.json
# start server
$WEST_ROOT/bin/w_run --work-manager=zmq --n-workers=0 --zmq-mode=master --zmq-write-host-info=$SERVER_INFO --zmq-comm-mode=tcp &> west-$SLURM_JOBID.log &

for ((n=0; n<60; n++)); do
    if [ -e $SERVER_INFO ] ; then
        echo "== server info file $SERVER_INFO =="
        cat $SERVER_INFO
        break
    fi
    sleep 1
done

# exit if host info file doesn't appear in one minute
if ! [ -e $SERVER_INFO ] ; then
    echo 'server failed to start'
    exit 1
fi

for node in $(scontrol show hostname $SLURM_NODELIST); do
   ssh -o StrictHostKeyChecking=no $node $PWD/node.sh $SLURM_SUBMIT_DIR $SLURM_JOBID $node --verbose \
                     --work-manager=zmq --n-workers=8 --zmq-mode=node \
                     --zmq-read-host-info=$SERVER_INFO \
                     &> west-$JOB_ID-$machine.log &
done

wait
