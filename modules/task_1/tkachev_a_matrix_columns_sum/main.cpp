// Copyright 2020 Tkachev Alexey
#include <gtest-mpi-listener.hpp>
#include <gtest/gtest.h>
#include <iostream>
#include <vector>
#include "../../../../modules/task_1/tkachev_a_matrix_columns_sum/tkachev_a_sum_matrix_columns.h"

TEST(Tkachev_Task1_MatrixColsSum, MatrixSize_20_20) {
    int rank, processes_count;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &processes_count);

    int count_columns = 20, count_rows = 20;
    std::vector<int> test_matrix(count_columns * count_rows);

    if (rank == 0) {
        test_matrix = randomMatrix(count_rows, count_columns);
    }

    std::vector<int> parallel_sum_cols(count_columns);
    parallel_sum_cols = parallelMatrixColumnsSum(test_matrix, count_columns, count_rows);

    if (rank == 0) {
        std::vector<int> sequential_sum(count_columns);
        sequential_sum = matrixColumnsSum(test_matrix, count_columns, count_rows,
                                                1, 0, test_matrix.size());
        ASSERT_EQ(parallel_sum_cols, sequential_sum);
    }
}

TEST(Tkachev_Task1_MatrixColsSum, MatrixSize_1_1) {
    int rank, processes_count;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &processes_count);

    int count_columns = 1, count_rows = 1;
    std::vector<int> test_matrix(count_columns * count_rows);

    if (rank == 0) {
        test_matrix = randomMatrix(count_rows, count_columns);
    }

    std::vector<int> parallel_sum_cols(count_columns);
    parallel_sum_cols = parallelMatrixColumnsSum(test_matrix, count_columns, count_rows);

    if (rank == 0) {
        std::vector<int> sequential_sum(count_columns);
        sequential_sum = matrixColumnsSum(test_matrix, count_columns, count_rows,
                                                1, 0, test_matrix.size());
        ASSERT_EQ(parallel_sum_cols, sequential_sum);
    }
}

TEST(Tkachev_Task1_MatrixColsSum, MatrixSize_1_10) {
    int rank, processes_count;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &processes_count);

    const int count_rows = 1;
    const int count_columns = 10;

    std::vector<int> test_matrix(count_columns * count_rows);

    if (rank == 0) {
        test_matrix = randomMatrix(count_rows, count_columns);
    }

    std::vector<int> parallel_sum_cols(count_columns);
    parallel_sum_cols = parallelMatrixColumnsSum(test_matrix, count_columns, count_rows);
    if (rank == 0) {
        std::vector<int> sequential_sum(count_columns);
        sequential_sum = matrixColumnsSum(test_matrix, count_columns, count_rows,
                                                1, 0, test_matrix.size());
        ASSERT_EQ(parallel_sum_cols, sequential_sum);
    }
}

TEST(Tkachev_Task1_MatrixColsSum, MatrixSize_31_1) {
    int rank, processes_count;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &processes_count);

    const int count_rows = 31;
    const int count_columns = 1;
    std::vector<int> test_matrix(count_columns * count_rows);

    if (rank == 0) {
        test_matrix = randomMatrix(count_rows, count_columns);
    }

    std::vector<int> parallel_sum_cols(count_columns);
    parallel_sum_cols = parallelMatrixColumnsSum(test_matrix, count_columns, count_rows);

    if (rank == 0) {
        std::vector<int> sequential_sum(count_columns);
        sequential_sum = matrixColumnsSum(test_matrix, count_columns, count_rows,
                                                1, 0, test_matrix.size());
        ASSERT_EQ(parallel_sum_cols, sequential_sum);
    }
}

TEST(Tkachev_Task1_MatrixColsSum, MatrixSize_4_3) {
    int rank, processes_count;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &processes_count);

    const int count_rows = 4;
    const int count_columns = 3;

    std::vector<int> test_matrix(count_columns * count_rows);

    if (rank == 0) {
        test_matrix = randomMatrix(count_rows, count_columns);
    }

    std::vector<int> parallel_sum_cols(count_columns);
    parallel_sum_cols = parallelMatrixColumnsSum(test_matrix, count_columns, count_rows);

    if (rank == 0) {
        std::vector<int> sequential_sum(count_columns);
        sequential_sum = matrixColumnsSum(test_matrix, count_columns, count_rows,
                                                1, 0, test_matrix.size());
        ASSERT_EQ(parallel_sum_cols, sequential_sum);
    }
}

TEST(Tkachev_Task1_MatrixColsSum, MatrixSize_21_20) {
    int rank, processes_count;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &processes_count);

    const int count_rows = 21;
    const int count_columns = 20;

    std::vector<int> test_matrix(count_columns * count_rows);

    if (rank == 0) {
        test_matrix = randomMatrix(count_rows, count_columns);
    }

    std::vector<int> parallel_sum_cols(count_columns);
    parallel_sum_cols = parallelMatrixColumnsSum(test_matrix, count_columns, count_rows);

    if (rank == 0) {
        std::vector<int> sequential_sum(count_columns);
        sequential_sum = matrixColumnsSum(test_matrix, count_columns, count_rows,
                                                1, 0, test_matrix.size());
        ASSERT_EQ(parallel_sum_cols, sequential_sum);
    }
}

TEST(Tkachev_Task1_MatrixColsSum, MatrixSize_41_41) {
    int rank, processes_count;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &processes_count);

    const int count_rows = 41;
    const int count_columns = 41;
    std::vector<int> test_matrix(count_columns * count_rows);

    if (rank == 0) {
        test_matrix = randomMatrix(count_rows, count_columns);
    }

    std::vector<int> parallel_sum_cols(count_columns);
    parallel_sum_cols = parallelMatrixColumnsSum(test_matrix, count_columns, count_rows);

    if (rank == 0) {
        std::vector<int> sequential_sum(count_columns);
        sequential_sum = matrixColumnsSum(test_matrix, count_columns, count_rows,
                                                1, 0, test_matrix.size());
        ASSERT_EQ(parallel_sum_cols, sequential_sum);
    }
}


int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    MPI_Init(&argc, &argv);

    ::testing::AddGlobalTestEnvironment(new GTestMPIListener::MPIEnvironment);
    ::testing::TestEventListeners& listeners =
            ::testing::UnitTest::GetInstance()->listeners();

    listeners.Release(listeners.default_result_printer());
    listeners.Release(listeners.default_xml_generator());

    listeners.Append(new GTestMPIListener::MPIMinimalistPrinter);
    return RUN_ALL_TESTS();
}
