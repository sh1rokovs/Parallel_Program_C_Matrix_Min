// Copyright 2020 Romanuyk Sergey
#ifndef MODULES_TASK_3_ROMANUYK_ALGORITM_FOX_ALGORITM_FOX_H_
#define MODULES_TASK_3_ROMANUYK_ALGORITM_FOX_ALGORITM_FOX_H_

#include <vector>

std::vector<double> getRandomMatrix(int orderM);
void BlockMult(double* pAblock, double* pBblock, double* pCblock, int blockSize);
std::vector<double> seqMult(const std::vector<double>& matA, const std::vector<double>& matB, int size);
std::vector<double> foxsAlgorithm(const std::vector<double>& matA, const std::vector<double>& matB, int orderM);
bool compareMat(const std::vector<double>& matrixA, const std::vector<double>& matrixB);

#endif  // MODULES_TASK_3_ROMANUYK_ALGORITM_FOX_ALGORITM_FOX_H_
