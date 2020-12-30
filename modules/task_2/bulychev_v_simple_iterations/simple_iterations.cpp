// Copyright 2020 Bulychev Vladislav
#include <mpi.h>
#include <vector>
#include <cmath>
#include "../../../modules/task_2/bulychev_v_simple_iterations/simple_iterations.h"

std::vector<double> Simple_iterations_Method(std::vector<double> mat,
    std::vector<double> vec, double eps) {
    int s = vec.size();
    if (mat.size() != s * s) {
        throw "error";
    }

    for (int i = 0; i < s; i++) {
        if (mat[i * s + i] == 0) {
            throw "error";
        }
    }

    std::vector<double> result(vec.size());

    for (int i = 0; i < s; i++) {
        result[i] = vec[i] / mat[i * s + i];
    }

    std::vector<double> x(s);

    do {
        for (int i = 0; i < s; i++) {
            x[i] = vec[i] / mat[i * s + i];
            for (int j = 0; j < s; j++) {
                if (i == j) {
                    continue;
                } else {
                    x[i] -= mat[i * s + j] / mat[i * s + i] * result[j];
                }
            }
        }

        bool flag = true;
        for (int i = 0; i < s - 1; i++) {
            if (std::abs(x[i] - result[i]) > eps) {
                flag = false;
                break;
            }
        }

        for (int i = 0; i < s; i++) {
            result[i] = x[i];
        }

        if (flag) {
            break;
        }
    } while (1);

    return result;
}

std::vector<double> MPI_Simple_iterations_Method(std::vector<double> mat,
    std::vector<double> vec, double eps) {
    int s = vec.size();
    int x1 = 0;

    if (mat.size() != s * s) {
        throw "error";
    }

    for (int i = 0; i < s; i++) {
        if (mat[i * s + i] == 0) {
            throw "error";
        }
    }

    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int lenght = s / size;

    if (x1 < size) {
        return Simple_iterations_Method(mat, vec, eps);
    }

    std::vector<double> l_mat;
    l_mat.resize(size);
    std::vector<double> l_vec;
    l_vec.resize(size);
    std::vector<double> tmp;
    tmp.resize(s);
    std::vector<double> x;
    x.resize(s);

    std::vector<int> counts_vec(size);
    std::vector<int> counts_mat(size);
    std::vector<int> displs_vec(size);
    std::vector<int> displs_mat(size);
    for (int i = 0; i < size; i++) {
        counts_vec[i] = lenght + s % size;
        counts_mat[i] = (lenght + s % size) * s;
        displs_vec[i] = lenght * i;
        displs_mat[i] = s % size + lenght * i;
    }

    MPI_Scatterv(&mat[0], &counts_mat[0], &displs_mat[0], MPI_DOUBLE,
        &l_mat[0], counts_mat[rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatterv(&vec[0], &counts_vec[0], &displs_vec[0], MPI_DOUBLE,
        &l_vec[0], counts_vec[rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);

    for (int i = 0; i < s; i++) {
        x[i] = l_vec[i] / l_mat[i * s + i];
    }

    MPI_Bcast(&x[0], s, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    int flag;
    do {
        for (int i = 0; i < s; i++) {
            int t = i * s + i;
            x[i] = vec[i] / mat[t];
            for (int j = 0; j < s; j++) {
                if (i == j) {
                    continue;
                } else {
                    x[i] -= l_mat[t] / mat[t] * tmp[j];
                }
            }
        }

        int local_flag = 1;
        for (int i = 0; i < s - 1; i++) {
            if (std::abs(x[i] - tmp[i]) > eps) {
                local_flag = 0;
                break;
            }
        }

        for (int i = 0; i < s; i++) {
            x[i] = tmp[i];
        }

        MPI_Gatherv(&tmp[0], counts_vec[rank], MPI_DOUBLE, &x[0],
            &counts_vec[0], &displs_vec[0], MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(&x[0], s, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        MPI_Allreduce(&local_flag, &flag, 1, MPI_INT, MPI_PROD, MPI_COMM_WORLD);
        if (flag == 1) {
            break;
        }
    } while (1);

    return x;
}
