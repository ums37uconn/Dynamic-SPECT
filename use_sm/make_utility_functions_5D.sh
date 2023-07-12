#!/bin/bash
#
# gcc flags
CFLAGS = "-g -O3 -m32"
# 64-bit file flages
RWFLAGS = "-D_FILE_OFFSET_BITS=64"
# OpenMP flag
OMPFLAGS = "-fopenmp"
echo ""
echo "Making utility binaries for LBL recon..."
echo ""
echo "Make sure an ./bin/ exists in this folder before running the program..."
echo ""
echo "gcc -g -O3 -m32 -D_FILE_OFFSET_BITS=64 EndianThis.c -o ./bin/EndianThis"
gcc -g -O3 -m32 -D_FILE_OFFSET_BITS=64 EndianThis.c -o ./bin/EndianThis
echo ""
echo "gcc -g -O3 -m32 -D_FILE_OFFSET_BITS=64 ManMat.c -o ./bin/ManMat"
gcc -g -O3 -m32 -D_FILE_OFFSET_BITS=64 ManMat.c -o ./bin/ManMat
echo ""
echo "gcc -g -O2 -m32 -D_FILE_OFFSET_BITS=64 MMD.c -o ./bin/MMD"
gcc -g -O3 -m32 -D_FILE_OFFSET_BITS=64 MMD.c -o ./bin/MMD
echo "gcc -g -O3 -m32 -D_FILE_OFFSET_BITS=64 SSR.c -lm -o ./bin/SSR"
gcc $-g -O3 -m32 -D_FILE_OFFSET_BITS=64 SSR.c -lm -o ./bin/SSR

# enabling 64-bit reads/writes
echo "gcc -g -O3 -m32 -D_FILE_OFFSET_BITS=64 use_sm.c -o ./bin/use_sm"
gcc -g -O3 -m32 -D_FILE_OFFSET_BITS=64 use_sm_5D.c -o ./use_sm_5D

# OpenMP version
# g++ $CFLAGS $RWFLAGS $OMPFLAGS use_sm_omp.cpp -o ./bin/use_sm_omp
# g++ -g -O2 -m32 -D_FILE_OFFSET_BITS=64 -fopenmp use_sm_omp.cpp -o ./bin/use_sm_omp
# echo "g++ -g -O2 -m32 -D_FILE_OFFSET_BITS=64 -fopenmp use_sm_omp.cpp -o ./bin/use_sm_omp"
