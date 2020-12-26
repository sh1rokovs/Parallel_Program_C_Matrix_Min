// Copyright 2020 Kasyanychev Mikhail
#ifndef MODULES_TASK_3_KASYANYCHEV_M_QUICK_SORT_EVEN_ODD_BATCHER_MERGE_KASYANYCHEV_M_QUICK_SORT_EVEN_ODD_BATCHER_MERGE_H_
#define MODULES_TASK_3_KASYANYCHEV_M_QUICK_SORT_EVEN_ODD_BATCHER_MERGE_KASYANYCHEV_M_QUICK_SORT_EVEN_ODD_BATCHER_MERGE_H_

#include <mpi.h>
#include <vector>

void getComp(std::vector<int>, std::vector<int>);
void getOddEven(std::vector<int>);
void createNet(int);
void genArray(std::vector<int>*);
void BatcherSort(std::vector<int>*);

#endif  // MODULES_TASK_3_KASYANYCHEV_M_QUICK_SORT_EVEN_ODD_BATCHER_MERGE_KASYANYCHEV_M_QUICK_SORT_EVEN_ODD_BATCHER_MERGE_H_
