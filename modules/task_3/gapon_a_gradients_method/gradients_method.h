// Copyright 2020 Gapon Andrey
#ifndef MODULES_TASK_3_GAPON_A_GRADIENTS_METHOD_GRADIENTS_METHOD_H_
#define MODULES_TASK_3_GAPON_A_GRADIENTS_METHOD_GRADIENTS_METHOD_H_
#include <vector>

std::vector<double> generate_random_vector(int size, double min, double max);
std::vector<double> generate_random_matrix(int size, double min, double max);
std::vector<double> gradients_method(const std::vector<double>& A, const std::vector<double>& b);
std::vector<double> gradients_method_parallel(const std::vector<double>& A, const std::vector<double>& b);

#endif  // MODULES_TASK_3_GAPON_A_GRADIENTS_METHOD_GRADIENTS_METHOD_H_
