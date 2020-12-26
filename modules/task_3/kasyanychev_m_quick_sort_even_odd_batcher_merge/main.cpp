// Copyright 2020 Kasyanychev Mikhail
#include <gtest/gtest.h>
#include <gtest-mpi-listener.hpp>
#include <algorithm>
#include <vector>
#include "./kasyanychev_m_quick_sort_even_odd_batcher_merge.h"

TEST(Batcher_Sort_MPI, Test_1) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::vector<int> arr{ 1, 7, 6, 4, 3, 8, 0, 9, 4, 5 };
    uint32_t sizeIn = arr.size();
    BatcherSort(&arr);
    if (rank == 0) {
        EXPECT_EQ(sizeIn, arr.size());
    }
}

TEST(Batcher_Sort_MPI, Test_2) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::vector<int> arrP{ 1, 9, 5, 6, 2, 3, 0, 4, 8, 6 };
    std::vector<int> arrS{ 1, 9, 5, 6, 2, 3, 0, 4, 8, 6 };
    BatcherSort(&arrP);
    if (rank == 0) {
        std::sort(arrS.begin(), arrS.end());
        bool AreEq = true;
        for (uint32_t i = 0; i < arrS.size(); i++) {
            if (arrS[i] != arrP[i]) {
                AreEq = false;
                break;
            }
        }
        EXPECT_EQ(true, AreEq);
    }
}

TEST(Batcher_Sort_MPI, Test_3) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::vector<int> arrP{ 1, 1, 1, 2, 2, 2, 0, 0, 0, 7 };
    std::vector<int> arrS{ 1, 1, 1, 2, 2, 2, 0, 0, 0, 7 };
    BatcherSort(&arrP);

    if (rank == 0) {
        std::sort(arrS.begin(), arrS.end());
        bool AreEq = true;
        for (uint32_t i = 0; i < arrS.size(); i++) {
            if (arrS[i] != arrP[i]) {
                AreEq = false;
                break;
            }
        }
        EXPECT_EQ(true, AreEq);
    }
}

TEST(Batcher_Sort_MPI, Test_4) {
    int rank, size = 1000;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::vector<int> arrS(size);

    if (rank == 0) {
        genArray(&arrS);
    }
    std::vector<int> arrP(arrS);
    BatcherSort(&arrP);
    if (rank == 0) {
        std::sort(arrS.begin(), arrS.end());
        bool AreEq = true;
        for (uint32_t i = 0; i < arrS.size(); i++) {
            if (arrS[i] != arrP[i]) {
                AreEq = false;
                break;
            }
        }
        EXPECT_EQ(true, AreEq);
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
