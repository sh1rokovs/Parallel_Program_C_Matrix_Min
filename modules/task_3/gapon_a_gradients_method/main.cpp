// Copyright 2020 Gapon Andrey
#include <gtest-mpi-listener.hpp>
#include <gtest/gtest.h>
#include <vector>
#include "./gradients_method.h"

TEST(Parallel_MPI, Test_Matr_1x1) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    std::vector<double> A = { 2.0 }, b = { 7.0 }, result = { 3.5 };
    const int N = 1;

    auto parallel_result = gradients_method_parallel(A, b);

    if (rank == 0) {
        auto seq_result = gradients_method(A, b);

        for (int i = 0; i < N; i++) {
            EXPECT_NEAR(parallel_result[i], result[i], 0.001);
            EXPECT_NEAR(seq_result[i], result[i], 0.001);
        }
    }
}

TEST(Parallel_MPI, Test_Matr_2x2) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    std::vector<double> A = { 4.56, 8.50, 8.50, 8.63 }, b = { 7.2, 5.8 }, result = { -0.39018, 1.05638 };
    const int N = 2;

    auto parallel_result = gradients_method_parallel(A, b);

    if (rank == 0) {
        auto seq_result = gradients_method(A, b);

        for (int i = 0; i < N; i++) {
            EXPECT_NEAR(parallel_result[i], result[i], 0.001);
            EXPECT_NEAR(seq_result[i], result[i], 0.001);
        }
    }
}

TEST(Parallel_MPI, Test_Matr_5x5) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    std::vector<double> A, b;
    const int N = 5;

    if (rank == 0) {
        A = generate_random_matrix(N, -10, 10);
        b = generate_random_vector(N, -10, 10);
    }

    auto parallel_result = gradients_method_parallel(A, b);

    if (rank == 0) {
        auto seq_result = gradients_method(A, b);

        for (int i = 0; i < N; i++) {
            EXPECT_NEAR(parallel_result[i], seq_result[i], 0.001);
        }
    }
}

TEST(Parallel_MPI, Test_Matr_15x15) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    std::vector<double> A, b;
    const int N = 15;

    if (rank == 0) {
        A = generate_random_matrix(N, -23, 17);
        b = generate_random_vector(N, -13, 53);
    }

    auto parallel_result = gradients_method_parallel(A, b);

    if (rank == 0) {
        auto seq_result = gradients_method(A, b);

        for (int i = 0; i < N; i++) {
            EXPECT_NEAR(parallel_result[i], seq_result[i], 0.001);
        }
    }
}

TEST(Parallel_MPI, Test_Matr_73x73) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    std::vector<double> A, b;
    const int N = 73;

    if (rank == 0) {
        A = generate_random_matrix(N, -108, 275);
        b = generate_random_vector(N, -94, 632);
    }

    auto parallel_result = gradients_method_parallel(A, b);

    if (rank == 0) {
        auto seq_result = gradients_method(A, b);

        for (int i = 0; i < N; i++) {
            EXPECT_NEAR(parallel_result[i], seq_result[i], 0.001);
        }
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    MPI_Init(&argc, &argv);

    ::testing::AddGlobalTestEnvironment(new GTestMPIListener::MPIEnvironment);
    ::testing::TestEventListeners& listeners =::testing::UnitTest::GetInstance()->listeners();

    listeners.Release(listeners.default_result_printer());
    listeners.Release(listeners.default_xml_generator());

    listeners.Append(new GTestMPIListener::MPIMinimalistPrinter);
    return RUN_ALL_TESTS();
}
