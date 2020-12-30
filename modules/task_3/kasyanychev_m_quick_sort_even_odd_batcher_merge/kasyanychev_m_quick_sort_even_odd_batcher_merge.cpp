// Copyright 2020 Kasyanychev Mikhail
#include "../../../modules/task_3/kasyanychev_m_quick_sort_even_odd_batcher_merge/kasyanychev_m_quick_sort_even_odd_batcher_merge.h"
#include <mpi.h>
#include <algorithm>
#include <ctime>
#include <iostream>
#include <random>
#include <utility>
#include <vector>

std::vector<std::pair<int, int>> comparators;

void genArray(std::vector<int>* array) {
    if (array->size() < 1) return;

    std::mt19937 gen;
    gen.seed(static_cast<unsigned int>(time(0)));

    for (uint32_t i = 0; i < array->size(); i++) {
        (*array)[i] = gen() % 64001 - 32000;
    }
    return;
}

void getComp(std::vector<int> left, std::vector<int> right) {
    int ressize = static_cast<int>(left.size()) + static_cast<int>(right.size());
    if (ressize == 1) return;
    if (ressize == 2) {
        std::pair<int, int> tmp{ left[0], right[0] };
        comparators.push_back(tmp);
        return;
    }
    std::vector<int> lvec_odd, rvec_odd, lvec_even, rvec_even, vecres(ressize);
    for (int i = 0; i < static_cast<int>(left.size()); i++) {
        if (i % 2)
            lvec_even.push_back(left[i]);
        else
            lvec_odd.push_back(left[i]);
    }
    for (int i = 0; i < static_cast<int>(right.size()); i++) {
        if (i % 2)
            rvec_even.push_back(right[i]);
        else
            rvec_odd.push_back(right[i]);
    }
    getComp(lvec_odd, rvec_odd);
    getComp(lvec_even, rvec_even);
    std::copy(left.begin(), left.end(), vecres.begin());
    std::copy(right.begin(), right.end(), vecres.begin() + left.size());
    for (int i = 1; i < static_cast<int>(vecres.size()) - 1; i += 2) {
        std::pair<int, int> tmp{ vecres[i], vecres[i + 1] };
        comparators.push_back(tmp);
    }
}

void getOddEven(std::vector<int> countP) {
    if (countP.size() < 2) return;
    std::vector<int> proc_left(countP.size() / 2);
    std::vector<int> proc_right(countP.size() / 2 + countP.size() % 2);
    std::copy(countP.begin(), countP.begin() + proc_left.size(), proc_left.begin());
    std::copy(countP.begin() + proc_left.size(), countP.end(), proc_right.begin());
    getOddEven(proc_left);
    getOddEven(proc_right);
    getComp(proc_left, proc_right);
}

void createNet(int size) {
    std::vector<int> proc(size);
    for (uint32_t i = 0; i < proc.size(); i++) {
        proc[i] = i;
    }
    getOddEven(proc);
}

void BatcherSort(std::vector<int>* res) {
    MPI_Status status;
    int rank, size;
    int _size = res->size();
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    createNet(size);
    int sizeNew = _size + ((_size % size) ? (size - (_size % size)) : 0);
    int elems_per_proc_size = sizeNew / size;
    for (int i = _size; i < sizeNew; i++) {
        res->push_back(INT16_MIN);
    }
    std::vector<int> elems_res(elems_per_proc_size);
    std::vector<int> elems_cur(elems_per_proc_size);
    std::vector<int> elems_tmp(elems_per_proc_size);
    MPI_Scatter(&(*res)[0], elems_per_proc_size, MPI_INT,
        &elems_res[0], elems_per_proc_size, MPI_INT, 0, MPI_COMM_WORLD);
    std::sort(elems_res.begin(), elems_res.end());
    for (uint32_t i = 0; i < comparators.size(); i++) {
        std::pair<int, int> comparator = comparators[i];
        if (rank == comparator.first) {
            MPI_Send(&elems_res[0], elems_per_proc_size, MPI_INT,
                comparator.second, 0, MPI_COMM_WORLD);
            MPI_Recv(&elems_cur[0], elems_per_proc_size, MPI_INT,
                comparator.second, 0, MPI_COMM_WORLD, &status);
            for (int resInd = 0, curInd = 0, tmpInd = 0;
                tmpInd < elems_per_proc_size; tmpInd++) {
                int res = elems_res[resInd];
                int cur = elems_cur[curInd];
                if (res < cur) {
                    elems_tmp[tmpInd] = res;
                    resInd++;
                } else {
                    elems_tmp[tmpInd] = cur;
                    curInd++;
                }
            }
            elems_res.swap(elems_tmp);
        } else if (rank == comparator.second) {
            MPI_Recv(&elems_cur[0], elems_per_proc_size, MPI_INT,
                comparator.first, 0, MPI_COMM_WORLD, &status);
            MPI_Send(&elems_res[0], elems_per_proc_size, MPI_INT,
                comparator.first, 0, MPI_COMM_WORLD);
            int start = elems_per_proc_size - 1;
            for (int resInd = start, curInd = start, tmpInd = start;
                tmpInd >= 0; tmpInd--) {
                int res = elems_res[resInd];
                int cur = elems_cur[curInd];
                if (res > cur) {
                    elems_tmp[tmpInd] = res;
                    resInd--;
                } else {
                    elems_tmp[tmpInd] = cur;
                    curInd--;
                }
            }
            elems_res.swap(elems_tmp);
        }
    }
    MPI_Gather(&elems_res[0], elems_per_proc_size, MPI_INT,
        &(*res)[0], elems_per_proc_size, MPI_INT, 0, MPI_COMM_WORLD);
    int diffElem = sizeNew - _size;
    if (rank == 0 && diffElem) {
        res->erase(res->begin(), res->begin() + diffElem);
    }
    return;
}
