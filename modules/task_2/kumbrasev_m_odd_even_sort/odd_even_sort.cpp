// Copyright 2020 Kumbrasev Mark
#include <mpi.h>
#include <iostream>
#include <cstdlib>
#include <random>
#include "../../../modules/task_2/kumbrasev_m_odd_even_sort/odd_even_sort.h"

long int gen_numbers(long int range) {
    std::random_device random;
    std::mt19937 generation(random());
    std::uniform_int_distribution<long int> uid(0, range);
    return uid(generation);
}

long int * gen_array(long int size) {
    long int * temp_arr = new long int[size];
    for (int i = 0; i < size; i++) {
        temp_arr[i] = gen_numbers(size);
    }
    return temp_arr;
}

void bubbleSort(long int *arr, long int size) {
     for (long int i = 0; i < size - 1; i++) {
        for (long int j = 0; j < size - i - 1; j++) {
            if (arr[j] > arr[j + 1]) {
                long int temp = arr[j];
                arr[j] = arr[j + 1];
                arr[j + 1] = temp;
            }
        }
    }
}

void odd_even_sort(long int * arr, long int size) {
    long int phase, itr;

    for (phase = 0; phase < size; phase++) {
        if (phase % 2 == 0) {
            for (itr = 1; itr < size; itr += 2) {
                if (arr[itr - 1] > arr[itr]) {
                    std::swap(arr[itr - 1], arr[itr]);
                }
            }
        }
        else {
            for (itr = 1; itr < size - 1; itr += 2) {
                if (arr[itr] > arr[itr + 1]) {
                    std::swap(arr[itr], arr[itr + 1]);
                }
            }
        }
    }
}

void PHASE(long int SEND_RANK, long int RCV_RANK, long int * arr, int size, MPI_Comm COMM) {
    int current_rank;
    MPI_Comm_rank(COMM, &current_rank);

    long int * temp_arr = new long int[size];

    long int * aux_arr = new long int[size * 2];

    if (current_rank == SEND_RANK) {
        MPI_Send(arr, size, MPI_LONG, RCV_RANK, 0, COMM);
        MPI_Recv(arr, size, MPI_LONG, RCV_RANK, 1, COMM, MPI_STATUS_IGNORE);
    }
    else {
        MPI_Recv(temp_arr, size, MPI_LONG, SEND_RANK, 0, COMM, MPI_STATUS_IGNORE);

        long int * first = &aux_arr[0];
        long int * last = &aux_arr[size * 2];

        long int * runner_1 = &arr[0];
        long int * runner_2 = &temp_arr[0];

        while (first != last) {
            if (runner_1 == &arr[size]) {
                while (runner_2 != &temp_arr[size]) {
                    *first++ = *runner_2++;
                }
            }
            else if (runner_2 == &temp_arr[size]) {
                while (runner_1 != &arr[size]) {
                    *first++ = *runner_1++;
                }
            }
            else if (*runner_1 < *runner_2) {
                *first++ = *runner_1++;
            }
            else {
                *first++ = *runner_2++;
            }
        }

        int itr = size;
        for (int i = 0; i < size; i++) {
            temp_arr[i] = aux_arr[i];
            arr[i] = aux_arr[itr];
            itr++;
        }

        delete[] aux_arr;

        MPI_Send(temp_arr, size, MPI_LONG, SEND_RANK, 1, COMM);

        delete[] temp_arr;
    }
}
