// Copyright 2020 Gapon Andrey
#ifndef MODULES_TASK_2_GAPON_A_GATHER_GATHER_H_
#define MODULES_TASK_2_GAPON_A_GATHER_GATHER_H_

#include <vector>

int Gather(void* sendbuf, int sendcount, MPI_Datatype sendtype,
    void* recbuf, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm);
std::vector<int> ArrInt(int num);
std::vector<float> ArrFloat(int num);
std::vector<double> ArrDouble(int num);

#endif  // MODULES_TASK_2_GAPON_A_GATHER_GATHER_H_
