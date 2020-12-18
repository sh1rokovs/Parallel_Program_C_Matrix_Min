// Copyright 2020 Romanuyk Sergey
#include <mpi.h>
#include <vector>
#include <random>
#include <ctime>
#include <numeric>
#include <iostream>
#include <algorithm>
#include <cstring>
#include <fstream>
#include <limits>
#include"../../../modules/task_3/romanuyk_algoritm_fox/algoritm_fox.h"

std::vector<double> genMatrix(int n) {
    std::mt19937 gen;
    gen.seed(static_cast<unsigned int>(time(0)));
    std::vector<double> arr(n*n);
    for (int i = 0; i < n * n; ++i) {
        arr[i] = gen() % 10 + static_cast<float>(gen() % 20) / 10;
    }
    return arr;
}

std::vector<double> SequentinalMultiMatrix(std::vector<double> A, std::vector<double> B, int n) {
    std::vector<double> res(n * n, 0);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            for (int k = 0; k < n; k++) {
                res[i * n + j] += A[i * n + k] * B[k * n + j];
            }
        }
    }
    return res;
}

void createGrid(Grid* grid, const int& procrank) {
    std::vector<int> GridCoords(2, 0);
    std::vector<int> DimSize(2, 0);
    std::vector<int> Periodic(2, 0);
    std::vector<int> Subdims(2, 0);
    DimSize[0] = grid->GridSize;
    DimSize[1] = grid->GridSize;
    Periodic[0] = 0;
    Periodic[1] = 0;

    MPI_Cart_create(MPI_COMM_WORLD, 2, DimSize.data(), Periodic.data(), 0, &grid->GridComm);

    MPI_Cart_coords(grid->GridComm, procrank, 2, GridCoords.data());
    grid->row = GridCoords[0];
    grid->col = GridCoords[1];

    Subdims[0] = 0;
    Subdims[1] = 1;
    MPI_Cart_sub(grid->GridComm, Subdims.data(), &grid->RowComm);

    Subdims[0] = 1;
    Subdims[1] = 0;
    MPI_Cart_sub(grid->GridComm, Subdims.data(), &grid->ColComm);
}

void MultiplyMatrixforParallel(const std::vector<double> A, const std::vector<double> matrixB,
    double* matrixC, int size) {
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            double temp = 0;
            for (int k = 0; k < size; k++)
                temp += A[i*size + k] * matrixB[k*size + j];
            matrixC[i*size + j] += temp;
        }
    }
}

std::vector<double> MultiplyMatrixParallel(const std::vector<double> A, const std::vector<double> B, int size) {
    int procnum, procrank;
    MPI_Comm_size(MPI_COMM_WORLD, &procnum);
    MPI_Comm_rank(MPI_COMM_WORLD, &procrank);
    MPI_Status status;

    Grid grid;
    grid.GridSize = static_cast<int>(sqrt(procnum));
    createGrid(&grid, procrank);

    int BlockSize;
    if (procrank == 0) {
        BlockSize = size / grid.GridSize;
    }
    MPI_Bcast(&BlockSize, 1, MPI_INT, 0, MPI_COMM_WORLD);

    MPI_Datatype matrixBlock;
    MPI_Type_vector(BlockSize, BlockSize, grid.GridSize * BlockSize, MPI_DOUBLE, &matrixBlock);
    MPI_Type_commit(&matrixBlock);

    std::vector<double> Ablock(BlockSize * BlockSize, 0);
    std::vector<double> Bblock(BlockSize * BlockSize, 0);

    if (procrank == 0) {
        for (int i = 0; i < BlockSize; i++) {
            for (int j = 0; j < BlockSize; j++) {
                Ablock[i * BlockSize + j] = A[i * size + j];
                Bblock[i * BlockSize + j] = B[i * size + j];
            }
        }
        for (int i = 1; i < procnum; i++) {
            int id = (i / grid.GridSize) * size * BlockSize + BlockSize * (i % grid.GridSize);
            MPI_Send(A.data() + id, 1, matrixBlock, i, 0, grid.GridComm);
            MPI_Send(B.data() + id, 1, matrixBlock, i, 1, grid.GridComm);
        }
    } else {
        MPI_Recv(Ablock.data(), BlockSize * BlockSize, MPI_DOUBLE, 0, 0, grid.GridComm, &status);
        MPI_Recv(Bblock.data(), BlockSize * BlockSize, MPI_DOUBLE, 0, 1, grid.GridComm, &status);
    }

    std::vector<double> Cblock(BlockSize * BlockSize, 0);

    for (int i = 0; i < grid.GridSize; i++) {
        std::vector<double> MatrixAblock(BlockSize * BlockSize, 0);
        int pivot = (grid.row + i) % grid.GridSize;
        if (grid.col == pivot) {
            MatrixAblock = Ablock;
        }
        MPI_Bcast(MatrixAblock.data(), BlockSize*BlockSize, MPI_DOUBLE, pivot, grid.RowComm);

        MultiplyMatrixforParallel(MatrixAblock, Bblock, Cblock.data(), BlockSize);

        int NextProc = grid.row + 1;
        if (grid.row == grid.GridSize - 1) NextProc = 0;
        int PrevProc = grid.row - 1;
        if (grid.row == 0) PrevProc = grid.GridSize - 1;
        MPI_Sendrecv_replace(Bblock.data(), BlockSize*BlockSize, MPI_DOUBLE, PrevProc, 0,
            NextProc, 0, grid.ColComm, &status);
    }

    std::vector<double> result(size*size);
    if (procrank == 0) {
        result.resize(size * size);
        for (int i = 0; i < BlockSize; i++) {
            for (int j = 0; j < BlockSize; j++) {
                result[i * size + j] = Cblock[i * BlockSize + j];
            }
        }
        for (int i = 1; i < procnum; i++) {
            MPI_Recv(result.data() + (i / grid.GridSize) * size * BlockSize + BlockSize * (i % grid.GridSize),
                1, matrixBlock, i, 3, MPI_COMM_WORLD, &status);
        }
    } else {
        MPI_Send(Cblock.data(), BlockSize * BlockSize, MPI_DOUBLE, 0, 3, MPI_COMM_WORLD);
    }
    MPI_Type_free(&matrixBlock);
    return result;
}

bool assertMatrix(const std::vector<double> A, const std::vector<double> B) {
    if (A.size() != B.size())
        throw "Different size";
    for (size_t i = 0; i < A.size(); i++) {
        if ((std::fabs(A[i] - B[i]) >= std::numeric_limits<double>::epsilon() * 1000000000.0))
            return false;
    }
    return true;
}
