// Copyright 2020 Ikromov Inom
#ifndef MODULES_TASK_2_IKROMOV_I_HYPERCUBE_HYPERCUBE_H_
#define MODULES_TASK_2_IKROMOV_I_HYPERCUBE_HYPERCUBE_H_

#include<mpi.h>

MPI_Comm getHypercube(int numDims = 0, int sizeNodeHyperCube = 0);
bool thisIsHypercube(MPI_Comm test_comm, int numDims = 0, int sizeNodeHyperCube = 0);
bool testHypercubeDataTransfer(MPI_Comm hypercube_comm, int numDims = 0, int sizeNodeHyperCube = 0);


#endif  // MODULES_TASK_2_IKROMOV_I_HYPERCUBE_HYPERCUBE_H_
