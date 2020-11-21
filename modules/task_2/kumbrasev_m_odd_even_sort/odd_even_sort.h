// Copyright 2020 Kumbrasev Mark
#ifndef MODULES_TASK_2_KUMBRASEV_M_ODD_EVEN_SORT_H_
#define MODULES_TASK_2_KUMBRASEV_M_ODD_EVEN_SORT_H_

#include <mpi.h>
#include <iostream>
#include <cstdlib>
#include <random>

long int gen_numbers(long int range);
long int * gen_array(long int size);

void bubbleSort(long int *arr, long int size);
void odd_even_sort(long int * arr, long int size);

void PHASE(long int SEND_RANK, long int RCV_RANK, long int * arr, int size, MPI_Comm COMM);

#endif  // MODULES_TASK_2_KUMBRASEV_M_ODD_EVEN_SORT_H_