// Copyright 2020 Gapon Andrey
#include <mpi.h>
#include <gtest-mpi-listener.hpp>
#include <gtest/gtest.h>
#include <iostream>
#include <vector>
#include "./gather.h"


TEST(Lab2, Test_Float) {
    int rank, size;
    int num = 800;
    int root = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    std::vector<float> vec = ArrFloat(num);
    std::vector<float> rbuf(size * num);
    std::vector<float> mybuf(size * num);
    Gather(vec.data(), num, MPI_FLOAT, mybuf.data(),
        num, MPI_FLOAT, root, MPI_COMM_WORLD);
    MPI_Gather(vec.data(), num, MPI_FLOAT, rbuf.data(),
        num, MPI_FLOAT, root, MPI_COMM_WORLD);
    if (rank == root) {
        ASSERT_EQ(mybuf, rbuf);
    }
}

TEST(Lab2, Test_Double) {
    int rank, size;
    int root = 0;
    int num = 800;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    std::vector<double> vec = ArrDouble(num);
    std::vector<double> rbuf(size* num);
    std::vector<double> mybuf(size* num);
    Gather(vec.data(), num, MPI_DOUBLE, mybuf.data(),
        num, MPI_DOUBLE, root, MPI_COMM_WORLD);
    MPI_Gather(vec.data(), num, MPI_DOUBLE, rbuf.data(),
        num, MPI_DOUBLE, root, MPI_COMM_WORLD);
    if (rank == root) {
        ASSERT_EQ(mybuf, rbuf);
    }
}


TEST(Lab2, Test_int_time) {
    int rank, size;
    int root = 0;
    int num = 1000;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    std::vector<int> vec = ArrInt(num);
    std::vector<int> rbuf(size * num);
    std::vector<int> mybuf(size * num);
    // auto start = std::chrono::steady_clock::now();
    Gather(vec.data(), num, MPI_INT, mybuf.data(),
        num, MPI_INT, root, MPI_COMM_WORLD);
    // auto end = std::chrono::steady_clock::now();
    // std::chrono::duration<double> my_timer = end - start;
    // auto start_mpi = std::chrono::steady_clock::now();
    MPI_Gather(vec.data(), num, MPI_INT, rbuf.data(),
        num, MPI_INT, root, MPI_COMM_WORLD);
    // auto end_mpi = std::chrono::steady_clock::now();
    // std::chrono::duration<double> mpi_timer = end_mpi - start_mpi;
    if (rank == root) {
        ASSERT_EQ(mybuf, rbuf);
    }
    if (rank == root) {
        // std::cout << "MPI guther runtime: "
        // << mpi_timer.count() << std::endl;
        // std::cout << "MPI guther runtime: "
        // << my_timer.count() << std::endl;
    }
}

TEST(Lab2, Test_double_time) {
    int rank, size;
    int root = 0;
    int num = 1000;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    std::vector<double> vec = ArrDouble(num);
    std::vector<double> rbuf(size * num);
    std::vector<double> mybuf(size * num);
    // auto start_2 = std::chrono::steady_clock::now();
    Gather(vec.data(), num, MPI_DOUBLE, mybuf.data(),
        num, MPI_DOUBLE, root, MPI_COMM_WORLD);
    // auto end_2 = std::chrono::steady_clock::now();
    // std::chrono::duration<double> my_timer = end_2 - start_2;
    // auto start_mpi_2 = std::chrono::steady_clock::now();
    MPI_Gather(vec.data(), num, MPI_DOUBLE, rbuf.data(),
        num, MPI_DOUBLE, root, MPI_COMM_WORLD);
    // auto end_mpi_2 = std::chrono::steady_clock::now();
    // std::chrono::duration<double> mpi_timer = end_mpi_2 - start_mpi_2;
    if (rank == root) {
        ASSERT_EQ(mybuf, rbuf);
    }
    if (rank == root) {
        // std::cout << "MPI guther runtime: "
        // << mpi_timer.count() << std::endl;
        // std::cout << "MPI guther runtime: "
        // << my_timer.count() << std::endl;
    }
}


int main(int argc, char* argv[]) {
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
