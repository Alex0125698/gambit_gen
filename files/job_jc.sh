#!/bin/bash 
#MSUB -r gambit_thdm            # Request name 
#MSUB -n 64                     # Total number of tasks to use 
#MSUB -c 1                      # Number of threads per task to use 
#MSUB -T 9999                   # Elapsed time limit in seconds (4 hours max)
#MSUB -o stdo_%I.txt              # Standard output. %I is the job id 
#MSUB -e stde_%I.txt              # Error output. %I is the job id
#MSUB -A ra5528                 # Project space
#MSUB -q rome                   # Partition
#MSUB -m scratch,work           # file systems


GAMBIT_DIR=${CCCSCRATCHDIR}/gambit
YAML_FILE=THDM_physical.yaml

source ~/scripts/loadModules.sh
set -x 
cd ${GAMBIT_DIR}
mkdir work
mkdir runs
rm -rf work/*
rm -rf runs/*
cd work

export OMP_NUM_THREADS=1
echo "Beginning scan on $(date '+%A %B %d %Y %r')"
start_time=$(date +%s)
timeout -s SIGTERM 8888 ccc_mprun ${GAMBIT_DIR}/gambit -rf ${GAMBIT_DIR}/yaml_files/${YAML_FILE}
cd ${GAMBIT_DIR}
echo "Finished scan on $(date '+%A %B %d %Y %r')"
end_time=$(date +%s)
echo "Job took $(( (end_time-start_time)/60 )) minutes."

start_time2=$(date +%s)
tar -c --use-compress-program=pigz -f scan.tar.gz ${GAMBIT_DIR}/runs/samples/scan.hdf5 && \
rm -rf ${GAMBIT_DIR}/runs/*
# rm -rf ${GAMBIT_DIR}/work/*
end_time2=$(date +%s)
echo "compression took $(( (end_time2-start_time2)/60 )) minutes."
