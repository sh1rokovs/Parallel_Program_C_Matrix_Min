// Copyright 2020 Gusarova Daria
#ifndef MODULES_TASK_3_GUSAROVA_D_ALG_STRONGIN_ALG_STRONGIN_H_
#define MODULES_TASK_3_GUSAROVA_D_ALG_STRONGIN_ALG_STRONGIN_H_

#include <mpi.h>
#include <functional>
#include <utility>

typedef std::pair<double, double> coords;
coords SequentalStrongin(const std::function<double(double)>& func, double a, double b);
coords ParallelStrongin(const std::function<double(double)>& func, double a, double b);


#endif  // MODULES_TASK_3_GUSAROVA_D_ALG_STRONGIN_ALG_STRONGIN_H_
