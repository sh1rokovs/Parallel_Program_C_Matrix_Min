// Copyright 2020 Ikromov Inom
#ifndef MODULES_TASK_2_IKROMOV_I_HYPERCUBE_HYPERCUBE_H_
#define MODULES_TASK_2_IKROMOV_I_HYPERCUBE_HYPERCUBE_H_

#include <mpi.h>

int getDimsCount(int proc_num);
MPI_Comm createHypercubeTopology(int dim_count);

#endif  // MODULES_TASK_2_IKROMOV_I_HYPERCUBE_HYPERCUBE_H_
