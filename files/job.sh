export OMP_NUM_THREADS=1
(mkdir work || true) && \
rm -rf work/* && \
cd work && \
mpiexec -np 16 ../gambit -rf ../yaml_files/THDM_physical.yaml > out.txt
