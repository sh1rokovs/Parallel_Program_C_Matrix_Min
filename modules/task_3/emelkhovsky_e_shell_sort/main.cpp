// Copyright 2020 Emelkhovsky Ekaterina
#include <gtest-mpi-listener.hpp>
#include <gtest/gtest.h>
#include <algorithm>
#include <vector>
#include "./shell_sort.h"

TEST(shell_sort, Test_1) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::vector<int> list = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    std::vector<int> result_list = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    std::vector<int> sortArray = shell_sort(list);
    if (rank == 0) {
        ASSERT_EQ(sortArray, result_list);
    }
}

TEST(shell_sort, Test_2) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::vector<int> list = {5, 6, 4, 3, 8, 2, 1, 10, 7, 9};
    std::vector<int> result_list = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    std::vector<int> sortArray = shell_sort(list);
    if (rank == 0) {
        ASSERT_EQ(sortArray, result_list);
    }
}

TEST(shell_sort, Test_3) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::vector<int> list = {-1, -2, -3, -4, -5};
    std::vector<int> result_list = {-5, -4, -3, -2, -1};
    std::vector<int> sortArray = shell_sort(list);
    if (rank == 0) {
        ASSERT_EQ(sortArray, result_list);
    }
}

TEST(shell_sort, Test_4) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::vector<int> list = {1, -10, 0, 123};
    std::vector<int> result_list = {-10, 0, 1, 123};
    std::vector<int> sortArray = shell_sort(list);
    if (rank == 0) {
        ASSERT_EQ(sortArray, result_list);
    }
}

TEST(shell_sort, Test_5) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::vector<int> list = {831};
    std::vector<int> result_list = {831};
    std::vector<int> sortArray = shell_sort(list);
    if (rank == 0) {
        ASSERT_EQ(sortArray, result_list);
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
