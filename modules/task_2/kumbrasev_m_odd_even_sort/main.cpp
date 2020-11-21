// Copyright 2020 Kumbrasev Mark
#include <gtest-mpi-listener.hpp>
#include <gtest/gtest.h>
#include <vector>
#include <algorithm>
#include "./odd_even_sort.h"

TEST(Parallel_Operations_MPI, Test_100) {
    int my_rank, comm_sz;
    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    long int size = 100;
    long int * arr = gen_array(size);
    int local_size = size / comm_sz;
    long int * local_arr = new long int[local_size];
    bubbleSort(arr, size);
    MPI_Scatter(arr, local_size, MPI_LONG, local_arr, local_size, MPI_LONG, 0, MPI_COMM_WORLD);
    delete[] arr;
    odd_even_sort(local_arr, local_size);
    for (int proc_itr = 1; proc_itr <= comm_sz; proc_itr++) {
        if ((my_rank + proc_itr) % 2 == 0) {
            if (my_rank < comm_sz - 1) {
                PHASE(my_rank, my_rank + 1, local_arr, local_size, MPI_COMM_WORLD);
            }
        }
        else if (my_rank > 0) {
            PHASE(my_rank - 1, my_rank, local_arr, local_size, MPI_COMM_WORLD);
        }
    }
    long int * final_arr = new long int[size];
    MPI_Gather(local_arr, local_size, MPI_LONG, final_arr, local_size, MPI_LONG, 0, MPI_COMM_WORLD);
    delete[] final_arr;
    if (my_rank == 0) {
    ASSERT_EQ(arr, final_arr);
    }
}

TEST(Parallel_Operations_MPI, Test_100) {
    int my_rank, comm_sz;
    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    long int size = 100;
    long int * arr = gen_array(size);
    int local_size = size / comm_sz;
    long int * local_arr = new long int[local_size];
    bubbleSort(arr, size);
    MPI_Scatter(arr, local_size, MPI_LONG, local_arr, local_size, MPI_LONG, 0, MPI_COMM_WORLD);
    delete[] arr;
    odd_even_sort(local_arr, local_size);
    for (int proc_itr = 1; proc_itr <= comm_sz; proc_itr++) {
        if ((my_rank + proc_itr) % 2 == 0) {
            if (my_rank < comm_sz - 1) {
                PHASE(my_rank, my_rank + 1, local_arr, local_size, MPI_COMM_WORLD);
            }
        }
        else if (my_rank > 0) {
            PHASE(my_rank - 1, my_rank, local_arr, local_size, MPI_COMM_WORLD);
        }
    }
    long int * final_arr = new long int[size];
    MPI_Gather(local_arr, local_size, MPI_LONG, final_arr, local_size, MPI_LONG, 0, MPI_COMM_WORLD);
    delete[] final_arr;
    if (my_rank == 0) {
    ASSERT_EQ(arr, final_arr);
    }
}

TEST(Parallel_Operations_MPI, Test_256) {
    int my_rank, comm_sz;
    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    long int size = 256;
    long int * arr = gen_array(size);
    int local_size = size / comm_sz;
    long int * local_arr = new long int[local_size];
    bubbleSort(arr, size);
    MPI_Scatter(arr, local_size, MPI_LONG, local_arr, local_size, MPI_LONG, 0, MPI_COMM_WORLD);
    delete[] arr;
    odd_even_sort(local_arr, local_size);
    for (int proc_itr = 1; proc_itr <= comm_sz; proc_itr++) {
        if ((my_rank + proc_itr) % 2 == 0) {
            if (my_rank < comm_sz - 1) {
                PHASE(my_rank, my_rank + 1, local_arr, local_size, MPI_COMM_WORLD);
            }
        }
        else if (my_rank > 0) {
            PHASE(my_rank - 1, my_rank, local_arr, local_size, MPI_COMM_WORLD);
        }
    }
    long int * final_arr = new long int[size];
    MPI_Gather(local_arr, local_size, MPI_LONG, final_arr, local_size, MPI_LONG, 0, MPI_COMM_WORLD);
    delete[] final_arr;
    if (my_rank == 0) {
    ASSERT_EQ(arr, final_arr);
    }
}

TEST(Parallel_Operations_MPI, Test_25) {
    int my_rank, comm_sz;
    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    long int size = 25;
    long int * arr = gen_array(size);
    int local_size = size / comm_sz;
    long int * local_arr = new long int[local_size];
    bubbleSort(arr, size);
    MPI_Scatter(arr, local_size, MPI_LONG, local_arr, local_size, MPI_LONG, 0, MPI_COMM_WORLD);
    delete[] arr;
    odd_even_sort(local_arr, local_size);
    for (int proc_itr = 1; proc_itr <= comm_sz; proc_itr++) {
        if ((my_rank + proc_itr) % 2 == 0) {
            if (my_rank < comm_sz - 1) {
                PHASE(my_rank, my_rank + 1, local_arr, local_size, MPI_COMM_WORLD);
            }
        }
        else if (my_rank > 0) {
            PHASE(my_rank - 1, my_rank, local_arr, local_size, MPI_COMM_WORLD);
        }
    }
    long int * final_arr = new long int[size];
    MPI_Gather(local_arr, local_size, MPI_LONG, final_arr, local_size, MPI_LONG, 0, MPI_COMM_WORLD);
    delete[] final_arr;
    if (my_rank == 0) {
    ASSERT_EQ(arr, final_arr);
    }
}

TEST(Parallel_Operations_MPI, Test_10) {
    int my_rank, comm_sz;
    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    long int size = 10;
    long int * arr = gen_array(size);
    int local_size = size / comm_sz;
    long int * local_arr = new long int[local_size];
    bubbleSort(arr, size);
    MPI_Scatter(arr, local_size, MPI_LONG, local_arr, local_size, MPI_LONG, 0, MPI_COMM_WORLD);
    delete[] arr;
    odd_even_sort(local_arr, local_size);
    for (int proc_itr = 1; proc_itr <= comm_sz; proc_itr++) {
        if ((my_rank + proc_itr) % 2 == 0) {
            if (my_rank < comm_sz - 1) {
                PHASE(my_rank, my_rank + 1, local_arr, local_size, MPI_COMM_WORLD);
            }
        }
        else if (my_rank > 0) {
            PHASE(my_rank - 1, my_rank, local_arr, local_size, MPI_COMM_WORLD);
        }
    }
    long int * final_arr = new long int[size];
    MPI_Gather(local_arr, local_size, MPI_LONG, final_arr, local_size, MPI_LONG, 0, MPI_COMM_WORLD);
    delete[] final_arr;
    if (my_rank == 0) {
    ASSERT_EQ(arr, final_arr);
    }
}

TEST(Parallel_Operations_MPI, Test_1250) {
    int my_rank, comm_sz;
    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    long int size = 1250;
    long int * arr = gen_array(size);
    int local_size = size / comm_sz;
    long int * local_arr = new long int[local_size];
    bubbleSort(arr, size);
    MPI_Scatter(arr, local_size, MPI_LONG, local_arr, local_size, MPI_LONG, 0, MPI_COMM_WORLD);
    delete[] arr;
    odd_even_sort(local_arr, local_size);
    for (int proc_itr = 1; proc_itr <= comm_sz; proc_itr++) {
        if ((my_rank + proc_itr) % 2 == 0) {
            if (my_rank < comm_sz - 1) {
                PHASE(my_rank, my_rank + 1, local_arr, local_size, MPI_COMM_WORLD);
            }
        }
        else if (my_rank > 0) {
            PHASE(my_rank - 1, my_rank, local_arr, local_size, MPI_COMM_WORLD);
        }
    }
    long int * final_arr = new long int[size];
    MPI_Gather(local_arr, local_size, MPI_LONG, final_arr, local_size, MPI_LONG, 0, MPI_COMM_WORLD);
    delete[] final_arr;
    if (my_rank == 0) {
    ASSERT_EQ(arr, final_arr);
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
