clear && printf '\033[3J'
echo "creating a release build of gambit..." && \
cd .. && \
(mkdir build || true ) && \
cd build && \
(make nuke-all && make clean || true ) && \
rm -rf ./* && \

# --- updated for JC below ---

export GAMBIT_DIR=${CCCSCRATCHDIR}/gambit
export BACKENDS_DIR=${GAMBIT_DIR}/Backends
export COPYME=${CCCSCRATCHDIR}/copyme

# copy pybind11 so that cmake won't download it
(mkdir -p ${GAMBIT_DIR}/contrib/pybind11/ || true) && \
cp -r ${COPYME}/contrib/pybind11/* ${GAMBIT_DIR}/contrib/pybind11

# cmake it!
cd ${GAMBIT_DIR}/build && \
cmake -D CMAKE_CXX_COMPILER=g++ -D CMAKE_C_COMPILER=gcc -D CMAKE_Fortran_COMPILER=gfortran -D CMAKE_BUILD_TYPE=Release_O3 -D WITH_MPI=ON -D BUILD_FS_MODELS="THDM_I;THDM_II;THDM_LS;THDM_flipped" -D WITH_ROOT=ON -D WITH_RESTFRAMES=ON -D WITH_HEPMC=ON -D EIGEN3_INCLUDE_DIR=/ccc/products/eigen-3.2.8/default/include/ .. && \
echo "\n\n\n\n"

# start the build process for each scanner / backend (it will fail)
timeout 2 make diver
touch diver_1.0.4-prefix/src/diver_1.0.4-stamp/diver_1.0.4-download
timeout 2 make THDMC
timeout 2 make higgsbounds
timeout 2 make heplikedata
timeout 2 make superiso
echo "\n\n\n\n"

# fudge the download stamp
touch THDMC_1.8.0-prefix/src/THDMC_1.8.0-stamp/THDMC_1.8.0-download
touch higgsbounds_5.8.0-prefix/src/higgsbounds_5.8.0-stamp/higgsbounds_5.8.0-download
touch heplikedata_1.0-prefix/src/heplikedata_1.0-stamp/heplikedata_1.0-download
touch superiso_4.1-prefix/src/superiso_4.1-stamp/superiso_4.1-download

# delete any old files
rm -rf ${GAMBIT_DIR}/ScannerBit/installed/diver/*
rm -rf ${BACKENDS_DIR}/scripts/BOSS/castxml/*
rm -rf ${BACKENDS_DIR}/installed/THDMC/1.8.0/*
rm -rf ${BACKENDS_DIR}/installed/higgsbounds/5.8.0/*
rm -rf ${BACKENDS_DIR}/installed/higgssignals/2.5.0/*
rm -rf ${BACKENDS_DIR}/installed/heplike/1.0/*
rm -rf ${BACKENDS_DIR}/installed/heplikedata/1.0/*
rm -rf ${BACKENDS_DIR}/installed/superiso/4.1/*

# copy the files across
(mkdir -p ${GAMBIT_DIR}/ScannerBit/installed/diver/ || true) && \
cp -r ${COPYME}/ScannerBit/installed/diver/* ${GAMBIT_DIR}/ScannerBit/installed/diver
(mkdir -p ${BACKENDS_DIR}/scripts/BOSS/castxml/ || true) && \
cp -r ${COPYME}/Backends/scripts/BOSS/castxml/* ${BACKENDS_DIR}/scripts/BOSS/castxml
(mkdir -p ${BACKENDS_DIR}/installed/THDMC/1.8.0/ || true) && \
cp -r ${COPYME}/Backends/installed/THDMC/1.8.0/* ${BACKENDS_DIR}/installed/THDMC/1.8.0
(mkdir -p ${BACKENDS_DIR}/installed/higgsbounds/5.8.0/build || true) && \
(mkdir -p ${BACKENDS_DIR}/installed/higgsbounds/5.8.0/data/lep-chisq-master/csboutput_trans_binary/ || true) && \
cp -r ${COPYME}/Backends/installed/higgsbounds/5.8.0/* ${BACKENDS_DIR}/installed/higgsbounds/5.8.0
(mkdir -p ${BACKENDS_DIR}/installed/heplikedata/1.0 || true) && \
cp -r ${COPYME}/Backends/installed/heplikedata/1.0/* ${BACKENDS_DIR}/installed/heplikedata/1.0
(mkdir -p ${BACKENDS_DIR}/installed/superiso/4.1 || true) && \
cp -r ${COPYME}/Backends/installed/superiso/4.1/* ${BACKENDS_DIR}/installed/superiso/4.1

# redo the build step
make diver
make THDMC
make higgsbounds
make heplikedata
make superiso

# higgssignals/heplike moved to the end
timeout 10 make higgssignals
touch higgssignals_2.5.0-prefix/src/higgssignals_2.5.0-stamp/higgssignals_2.5.0-download
rm -rf ${BACKENDS_DIR}/installed/higgssignals/2.5.0/*
(mkdir -p ${BACKENDS_DIR}/installed/higgssignals/2.5.0/build || true) && \
cp -r ${COPYME}/Backends/installed/higgssignals/2.5.0/* ${BACKENDS_DIR}/installed/higgssignals/2.5.0
make higgssignals

timeout 2 make heplike
touch heplike_1.0-prefix/src/heplike_1.0-stamp/heplike_1.0-download
rm -rf ${BACKENDS_DIR}/installed/heplike/1.0/*
(mkdir -p ${BACKENDS_DIR}/installed/heplike/1.0 || true) && \
cp -r ${COPYME}/Backends/installed/heplike/1.0/* ${BACKENDS_DIR}/installed/heplike/1.0

cd ${BACKENDS_DIR}/installed/heplike/1.0 &&
cmake -DCMAKE_CXX_FLAGS="-I/ccc/products/boost-1.69.0/gcc--8.3.0__openmpi--4.0.1/default/include/ -I/ccc/products/yaml-cpp-0.6.2/system/default/include" .
cd ${GAMBIT_DIR}/build
make heplike

# build gambit
cd ${GAMBIT_DIR}/build
cmake ..
make -j64 gambit
