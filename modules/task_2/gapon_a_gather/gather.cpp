// Copyright 2020 Gapon Andrey
#include <mpi.h>
#include <ctime>
#include <cstring>
#include <algorithm>
#include <random>
#include <iostream>
#include <vector>
#include <string>
#include "../../../modules/task_2/gapon_a_gather/gather.h"

std::vector<int> ArrInt(int size) {
    std::mt19937 gen;
    gen.seed(static_cast<unsigned int>(time(0)));

    std::vector<int> array(size);
    for (auto itr_int = array.begin(); itr_int !=array.end(); itr_int++)
        *itr_int = (gen() % 200);

    return array;
}

std::vector<float> ArrFloat(int size) {
    std::mt19937 gen;
    gen.seed(static_cast<unsigned int>(time(0)));

    std::vector<float> array(size);
    for (auto itr_fl = array.begin(); itr_fl != array.end(); itr_fl++) {
        *itr_fl = ((gen() % 200) / 10.0);
    }

    return array;
}

std::vector<double> ArrDouble(int size) {
    std::mt19937 gen;
    gen.seed(static_cast<unsigned int>(time(0)));

    std::vector<double> array(size);
    for (auto itr_d = array.begin(); itr_d != array.end(); itr_d++) {
      *itr_d = ((gen() % 200) / 10.0);
    }

    return array;
}

int Gather(void *sendbuf, int sendcount, MPI_Datatype sendtype,
    void *recbuf, int recvcount, MPI_Datatype recvtype,
    int root, MPI_Comm comm) {
    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (root < 0 || root >= size)
        return MPI_ERR_ROOT;
    int send_type_size, rec_type_size;
    MPI_Type_size(sendtype, &send_type_size);
    MPI_Type_size(recvtype, &rec_type_size);
    if (rec_type_size == MPI_ERR_TYPE || send_type_size == MPI_ERR_TYPE)
        return MPI_ERR_TYPE;
    if (sendcount <= 0 || recvcount <= 0 || sendcount != recvcount)
        return MPI_ERR_COUNT;
    if (rank == root) {
        for (int j = 0; j < size; j++) {
            if (j != root) {
                MPI_Status status;
                MPI_Recv(reinterpret_cast<char*>(recbuf) + recvcount *
                    rec_type_size * j, recvcount,
                    recvtype, j, 0, comm, &status);
            } else {
                std::memcpy(reinterpret_cast<char*>(recbuf) +
                    recvcount * rec_type_size * root,
                    sendbuf, sendcount* send_type_size);
            }
        }
    } else {
        MPI_Send(sendbuf, sendcount, sendtype, root, 0, comm);
    }
    return MPI_SUCCESS;
}
