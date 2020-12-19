// Copyright 2020 Romanuyk Sergey
#ifndef MODULES_TASK_3_ROMANUYK_ALGORITM_FOX_ALGORITM_FOX_H_
#define MODULES_TASK_3_ROMANUYK_ALGORITM_FOX_ALGORITM_FOX_H_

#include <vector>

std::vector<double> genMatrix(int size);
std::vector<double> SequentinalMultiMatrix(std::vector<double>& A, std::vector<double>& B, int size);

bool assertMatrix(const std::vector<double>& A, const std::vector<double>& B);

void createGrid(int gridSize, int ProcRank, int* gridCoords);
void MultiplyMatrixforParallel(double* pAblock, double* pBblock, double* pCblock, int blockSize);

std::vector<double> MultiplyMatrixParallel(const std::vector<double>& matA, const std::vector<double>& matB, int size);

#endif  // MODULES_TASK_3_ROMANUYK_ALGORITM_FOX_ALGORITM_FOX_H_
