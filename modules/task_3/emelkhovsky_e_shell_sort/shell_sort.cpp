// Copyright 2020 Emelkhovsky Ekaterina

#include <mpi.h>
#include <algorithm>
#include <iostream>
#include <vector>
#include "../../../modules/task_3/emelkhovsky_e_shell_sort/shell_sort.h"

std::vector<int> compareResults(std::vector<int> list) {
    int len = list.size();

    for (int i = 1; i < len; i++) {
        int j = 1;
        while ((i - j >= 0) && (list[i - j + 1] < list[i - j])) {
            int tmp = list[i - j + 1];
            list[i - j + 1] = list[i - j];
            list[i - j] = tmp;
            j++;
        }
    }
    return list;
}

std::vector<int> shell_sort(std::vector<int> list) {
    if (list.empty() || (list.size() < 2))
        return list;

    int procNum;
    MPI_Comm_size(MPI_COMM_WORLD, &procNum);
    int procRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &procRank);
    MPI_Status status;

    int len = list.size();
    int part = len;
    int lenin, count_of_parts, count_of_processes, count, count_2;
    int tmp = 0;
    int countProcNum = 0;


    while (part != 0) {
        part = part / 2;
        count_of_parts = part / procNum;
        count_of_processes = part % procNum;

        if (part < procNum) {
            count = part;
        } else {
            count = procNum;
        }

        count_2 = count_of_parts;

        if (count_of_processes > 0) {
            count_2++;
        }

        tmp = 0;
        countProcNum = 0;
        lenin = len / part;
        std::vector<int> list_value(lenin);
        std::vector<int> localArrayForRoot(lenin);

        do {
            if (countProcNum + count_of_processes == part) {
                count = count_of_processes;
            }

            for (int proc = 0; proc < count; proc++) {
                if (procRank == 0) {
                    list_value.clear();

                    for (int i = 0; i < lenin; i++) {
                        list_value.push_back(list[proc + countProcNum + part * i]);
                    }

                    if (proc == 0) {
                        list_value = compareResults(list_value);
                        for (int i = 0; i < lenin; i++) {
                            list[countProcNum + part * i] = list_value[i];
                        }
                    } else {
                        MPI_Send(&list_value[0], lenin, MPI_INT, proc, tmp, MPI_COMM_WORLD);
                    }
                } else {
                    if (procRank == proc) {
                        MPI_Recv(&list_value[0], lenin, MPI_INT, 0, tmp, MPI_COMM_WORLD, &status);
                        list_value = compareResults(list_value);
                        MPI_Send(&list_value[0], lenin, MPI_INT, 0, proc + tmp, MPI_COMM_WORLD);
                    }
                }
            }

            countProcNum = countProcNum + procNum;
            tmp++;
        } while (countProcNum < part);
        if (part < procNum) {
            count = part;
        } else {
            count = procNum;
        }

        countProcNum = 0;
        tmp = 0;

        if (procRank == 0) {
            do {
                if (countProcNum + count_of_processes == part)
                    count = count_of_processes;
                for (int proc = 1; proc < count; proc++) {
                    MPI_Recv(&list_value[0], lenin, MPI_INT, proc, proc + tmp, MPI_COMM_WORLD, &status);
                    for (int i = 0; i < lenin; i++) {
                        list[proc + part * i + countProcNum] = list_value[i];
                    }
                }
                countProcNum = countProcNum + procNum;
                tmp++;
            } while (countProcNum < part);
        }
        if (lenin == len) {
            part = 0;
        }
    }
    return list;
}
