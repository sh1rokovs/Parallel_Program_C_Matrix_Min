// Copyright 2020 Gapon Andrey
#include <mpi.h>
#include <random>
#include <vector>
#include <cmath>
#include <iostream>
#include "../../../modules/task_3/gapon_a_gradients_method/gradients_method.h"

std::vector<double> generate_random_vector(int size, double min, double max) {
    std::random_device random;
    std::mt19937 generator;
    generator.seed(random());
    std::uniform_real_distribution<> distr(min, max);

    std::vector<double> result(size);
    for (auto& el : result) {
        el = std::trunc(distr(generator));
    }

    return result;
}

std::vector<double> generate_random_matrix(int size, double min, double max) {
    std::random_device random;
    std::mt19937 generator;
    generator.seed(random());
    std::uniform_real_distribution<> distr(min, max);

    std::vector<double> result(size * size);
    std::vector<double> result2(size * size);

    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            result[i * size + j] = std::trunc(distr(generator));
        }
    }

    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            result2[i * size + j] = (result[i * size + j] = result[j * size + i]) / 2.0;
        }
    }

    return result2;
}

double scalar_product(const std::vector<double>& vec1, const std::vector<double>& vec2) {
    double result = 0.0;

    for (auto iter1 = vec1.cbegin(), iter2 = vec2.cbegin(); iter1 < vec1.end(); ++iter1, ++iter2) {
        result += *iter1 * *iter2;
    }

    return result;
}

std::vector<double> matr_vec_mult(const std::vector<double>& matr, const std::vector<double>& vec) {
    int size = static_cast<int>(vec.size());
    std::vector<double> result(size);

    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++)
            result[i] += matr[i * size + j] * vec[j];
    }

    return result;
}

std::vector<double> vec_matr_mult(const std::vector<double>& vec, const std::vector<double>& matr) {
    int size = static_cast<int>(vec.size());
    std::vector<double> result(size);

    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++)
            result[i] += vec[j] * matr[j * size + i];
    }

    return result;
}

std::vector<double> vec_diff(const std::vector<double>& vec1, const std::vector<double>& vec2) {
    int size = static_cast<int>(vec1.size());
    std::vector<double> result(size);

    for (int i = 0; i < size; i++) {
        result[i] = vec1[i] - vec2[i];
    }

    return result;
}

std::vector<double> vec_sum(const std::vector<double>& vec1, const std::vector<double>& vec2) {
    int size = static_cast<int>(vec1.size());
    std::vector<double> result(size);

    for (int i = 0; i < size; i++) {
        result[i] = vec1[i] + vec2[i];
    }

    return result;
}

std::vector<double> scalar_vec_mult(const double scalar, const std::vector<double>& vec) {
    int size = static_cast<int>(vec.size());
    std::vector<double> result(size);

    for (int i = 0; i < size; i++) {
        result[i] = scalar * vec[i];
    }

    return result;
}

std::vector<double> gradients_method(const std::vector<double>& A, const std::vector<double>& b) {
    int size = static_cast<double>(b.size());

    std::vector<double> x(size), w(size);

    auto z = b, r = b;
    double s = scalar_product(r, r);

    for (int i = 0; i < size; ++i) {
        w = matr_vec_mult(A, z);

        double m = scalar_product(w, z);
        double alpha = s / m;
        for (int j = 0; j < size; ++j) {
            x[j] += alpha * z[j];
            r[j] -= alpha * w[j];
        }

        double m2 = scalar_product(r, r);
        double beta = m2 / s;
        s = m2;
        for (int j = 0; j < size; ++j) {
            z[j] = beta * z[j] + r[j];
        }
    }

    return x;
}

std::vector<double> gradients_method_parallel(const std::vector<double>& A, const std::vector<double>& b) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int matr_size = static_cast<int>(b.size());

    MPI_Bcast(&matr_size, 1, MPI_INT, 0, MPI_COMM_WORLD);

    int delta = matr_size / size;
    int remainder = matr_size % size;

    std::vector<int> displaces(size), send_count(size), rows_count(size), displaces_matr(size);
    int displaces_count = 0;

    for (int i = 0; i < size; i++) {
        if (i < remainder)
            send_count[i] = delta + 1;
        else
            send_count[i] = delta;

        rows_count[i] = send_count[i] * matr_size;

        displaces[i] = displaces_count;
        displaces_matr[i] = displaces[i] * matr_size;
        displaces_count += send_count[i];
    }

    std::vector<double> localA(rows_count[rank]);
    std::vector<double> localb(send_count[rank]), localB(b), res(send_count[rank]);

    MPI_Scatterv(A.data(), rows_count.data(), displaces_matr.data(), MPI_DOUBLE,
        localA.data(), rows_count[rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatterv(b.data(), send_count.data(), displaces.data(), MPI_DOUBLE,
        localb.data(), send_count[rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
    if (rank != 0) {
        localB.resize(matr_size);
    }
    MPI_Bcast(localB.data(), matr_size, MPI_INT, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        res = gradients_method(A, localB);
    } else {
        res = gradients_method(localA, localb);
    }
    std::vector<double> result;

    if (rank == 0)
        result.resize(matr_size);

    MPI_Gatherv(res.data(), send_count[rank], MPI_DOUBLE,
        result.data(), send_count.data(), displaces.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD);

    return res;
}
