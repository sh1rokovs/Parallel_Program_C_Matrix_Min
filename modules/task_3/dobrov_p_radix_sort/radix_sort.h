// Copyright 2020 Dobrov Pavel

#ifndef MODULES_TASK_3_DOBROV_P_RADIX_SORT_RADIX_SORT_H_
#define MODULES_TASK_3_DOBROV_P_RADIX_SORT_RADIX_SORT_H_

#include <vector>

typedef struct {
    int rank1;
    int rank2;
} pair;

std::vector<double> getRandVec(int size);

std::vector<double> sortingByCounting(std::vector<double> vec1,
    std::vector<double> vec2, int byte);
std::vector<double> radixSort(std::vector<double> vec);

void batcher(int countOfProc);
void buildNetwork(std::vector<int> prcsVec);
void buildConnection(std::vector<int> upPrcsVec,
    std::vector<int> downPrcsVec);

std::vector<double> parOddEvenMerge(std::vector<double> globalVec);

#endif  // MODULES_TASK_3_DOBROV_P_RADIX_SORT_RADIX_SORT_H_
