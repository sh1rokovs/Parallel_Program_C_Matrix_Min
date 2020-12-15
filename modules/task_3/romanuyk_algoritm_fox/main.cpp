// Copyright 2020 Romanuyk Sergey
#include <gtest-mpi-listener.hpp>
#include <gtest/gtest.h>
#include <vector>
#include <ctime>
#include"../../../modules/task_3/romanuyk_algoritm_fox/algoritm_fox.h"

TEST(Parallel_Operations_MPI, Test1) {
    int procrank;
    MPI_Comm_rank(MPI_COMM_WORLD, &procrank);
    std::vector<double> a;
    std::vector<double> b;
    int n = 4;
    if (procrank == 0) {
        a = genMatrix(n);
        b = genMatrix(n);
    }
    double t_b1, t_e1, t_b2, t_e2;
    if (procrank == 0) {
        t_b1 = MPI_Wtime();
    }
    std::vector<double> res1 = MultiplyMatrixParallel(a, b, n);
    if (procrank == 0) {
        t_e1 = MPI_Wtime();
        std::cout << "Parallel time: " << t_e1 - t_b1 << std::endl;
    }
    if (procrank == 0) {
        t_b2 = MPI_Wtime();
        std::vector<double> res2 = SequentinalMultiMatrix(a, b, n);
        t_e2 = MPI_Wtime();
        std::cout << "Sequentional time: " << t_e2 - t_b2 << std::endl;
        ASSERT_TRUE(assertMatrix(res1, res2));
    }
}

TEST(Parallel_Operations_MPI, Test2) {
    int procrank;
    MPI_Comm_rank(MPI_COMM_WORLD, &procrank);
    std::vector<double> a;
    std::vector<double> b;
    int n = 16;
    if (procrank == 0) {
        a = genMatrix(n);
        b = genMatrix(n);
    }
    double t_b1, t_e1, t_b2, t_e2;
    if (procrank == 0) {
        t_b1 = MPI_Wtime();
    }
    std::vector<double> res1 = MultiplyMatrixParallel(a, b, n);
    if (procrank == 0) {
        t_e1 = MPI_Wtime();
        std::cout << "Parallel time: " << t_e1 - t_b1 << std::endl;
    }
    if (procrank == 0) {
        t_b2 = MPI_Wtime();
        std::vector<double> res2 = SequentinalMultiMatrix(a, b, n);
        t_e2 = MPI_Wtime();
        std::cout << "Sequentional time: " << t_e2 - t_b2 << std::endl;
        ASSERT_TRUE(assertMatrix(res1, res2));
    }
}

TEST(Parallel_Operations_MPI, Test3) {
    int procrank;
    MPI_Comm_rank(MPI_COMM_WORLD, &procrank);
    std::vector<double> a;
    std::vector<double> b;
    int n = 32;
    if (procrank == 0) {
        a = genMatrix(n);
        b = genMatrix(n);
    }
    double t_b1, t_e1, t_b2, t_e2;
    if (procrank == 0) {
        t_b1 = MPI_Wtime();
    }
    std::vector<double> res1 = MultiplyMatrixParallel(a, b, n);
    if (procrank == 0) {
        t_e1 = MPI_Wtime();
        std::cout << "Parallel time: " << t_e1 - t_b1 << std::endl;
    }
    if (procrank == 0) {
        t_b2 = MPI_Wtime();
        std::vector<double> res2 = SequentinalMultiMatrix(a, b, n);
        t_e2 = MPI_Wtime();
        std::cout << "Sequentional time: " << t_e2 - t_b2 << std::endl;
        ASSERT_TRUE(assertMatrix(res1, res2));
    }
}

TEST(Parallel_Operations_MPI, Test4) {
    int procrank;
    MPI_Comm_rank(MPI_COMM_WORLD, &procrank);
    std::vector<double> a;
    std::vector<double> b;
    int n = 100;
    if (procrank == 0) {
        a = genMatrix(n);
        b = genMatrix(n);
    }
    double t_b1, t_e1, t_b2, t_e2;
    if (procrank == 0) {
        t_b1 = MPI_Wtime();
    }
    std::vector<double> res1 = MultiplyMatrixParallel(a, b, n);
    if (procrank == 0) {
        t_e1 = MPI_Wtime();
        std::cout << "Parallel time: " << t_e1 - t_b1 << std::endl;
    }
    if (procrank == 0) {
        t_b2 = MPI_Wtime();
        std::vector<double> res2 = SequentinalMultiMatrix(a, b, n);
        t_e2 = MPI_Wtime();
        std::cout << "Sequentional time: " << t_e2 - t_b2 << std::endl;
        ASSERT_TRUE(assertMatrix(res1, res2));
    }
}

TEST(Parallel_Operations_MPI, Test5) {
    int procrank;
    MPI_Comm_rank(MPI_COMM_WORLD, &procrank);
    std::vector<double> a;
    std::vector<double> b;
    int n = 64;
    if (procrank == 0) {
        a = genMatrix(n);
        b = genMatrix(n);
    }
    double t_b1, t_e1, t_b2, t_e2;
    if (procrank == 0) {
        t_b1 = MPI_Wtime();
    }
    std::vector<double> res1 = MultiplyMatrixParallel(a, b, n);
    if (procrank == 0) {
        t_e1 = MPI_Wtime();
        std::cout << "Parallel time: " << t_e1 - t_b1 << std::endl;
    }
    if (procrank == 0) {
        t_b2 = MPI_Wtime();
        std::vector<double> res2 = SequentinalMultiMatrix(a, b, n);
        t_e2 = MPI_Wtime();
        std::cout << "Sequentional time: " << t_e2 - t_b2 << std::endl;
        ASSERT_TRUE(assertMatrix(res1, res2));
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
