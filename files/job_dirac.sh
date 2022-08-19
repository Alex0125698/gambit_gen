#!/bin/bash

#! NOTES
#! see hours using mybalance
#! submit using: sbatch submission_script
#! interactive job using: sintr -A DIRAC-DP192-CPU -p icelake -N1 -n76 -t 0:10:0 --qos=INTR
#! interactive job using: sintr -A DIRAC-DP192-CPU -p icelake -N1 -n38 -t 0:10:0 --qos=INTR
#! view jobs using: 
#! cancel using: scancel  id
#! cores/threads = 38/76 *2 (NB hyperthreading is actually disabled)
#! The --ntasks value refers to the number of tasks to be launched by SLURM only. This
#! usually equates to the number of MPI tasks launched. Reduce this from nodes*76 if
#! demanded by memory requirements, or if OMP_NUM_THREADS>1.

#!#############################################################

#SBATCH -J gambit_thdm
#SBATCH -A DIRAC-DP192-CPU
#SBATCH -p icelake
#SBATCH --nodes=1
#! How many (MPI) tasks will there be in total? (<= nodes*76)
#! The Ice Lake (icelake) nodes have 76 CPUs (cores) each and
#! 3380 MiB of memory per CPU.
#SBATCH --ntasks=76
#SBATCH --time=HHHHH:MMMMM:SSSSS
#SBATCH --mail-type=NONE
##SBATCH --no-requeue

numnodes=$SLURM_JOB_NUM_NODES # don't touch
numtasks=$SLURM_NTASKS # don't touch
mpi_tasks_per_node=$(echo "$SLURM_TASKS_PER_NODE" | sed -e  's/^\([0-9][0-9]*\).*$/\1/') # don't touch
. /etc/profile.d/modules.sh # don't touch

# set options below
source ~/scripts/loadModules.sh
export OMP_NUM_THREADS=1
application="../gambit"
options="-rf ../yaml_files/THDM_physical.yaml"

np=$[${numnodes}*${mpi_tasks_per_node}] # don't touch
export I_MPI_PIN_DOMAIN=omp:compact # don't touch
export I_MPI_PIN_ORDER=scatter # don't touch

#! Choose this for a MPI code (possibly using OpenMP) using Intel MPI.
# CMD="mpirun -ppn $mpi_tasks_per_node -np $np $application $options"

#! Choose this for a pure shared-memory OpenMP parallel program on a single node:
# CMD="$application $options"

#! Choose this for a MPI code (possibly using OpenMP) using OpenMPI:
CMD="mpirun -npernode $mpi_tasks_per_node -np $np $application $options"

# setup work dir below
cd GAMBIT_BASE_DIR
(mkdir work || true) && \
rm -rf work/* && \
cd work

# dont touch below
JOBID=$SLURM_JOB_ID

echo -e "JobID: $JOBID\n======"
echo "Time: `date`"
echo "Running on master node: `hostname`"
echo "Current directory: `pwd`"

#! Create a machine file:
if [ "$SLURM_JOB_NODELIST" ]; then
        export NODEFILE=`generate_pbs_nodefile`
        cat $NODEFILE | uniq > machine.file.$JOBID
        echo -e "\nNodes allocated:\n================"
        echo `cat machine.file.$JOBID | sed -e 's/\..*$//g'`
fi

echo -e "\nnumtasks=$numtasks, numnodes=$numnodes, mpi_tasks_per_node=$mpi_tasks_per_node (OMP_NUM_THREADS=$OMP_NUM_THREADS)"
echo -e "\nExecuting command:\n==================\n$CMD\n"

# eval $CMD 

echo -e "Beginning scan on $(date '+%A %B %d %Y %r')"
start_time=$(date +%s)
timeout -s SIGTERM 878787 mpirun -npernode $mpi_tasks_per_node -np $np $application $options >out.txt 2>err.txt
echo -e "Finished scan on $(date '+%A %B %d %Y %r')"
end_time=$(date +%s)
echo -e "Job took $(( (end_time-start_time)/60 )) minutes."

# start_time2=$(date +%s)
# cd ${GAMBIT_DIR}
# tar -c --use-compress-program=pigz -f scan.tar.gz ${GAMBIT_DIR}/runs/samples/scan.hdf5 && \
# rm -rf ${GAMBIT_DIR}/runs/*
# # rm -rf ${GAMBIT_DIR}/work/*
# end_time2=$(date +%s)
# echo "compression took $(( (end_time2-start_time2)/60 )) minutes."