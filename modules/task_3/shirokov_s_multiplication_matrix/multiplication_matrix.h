// Copyright 2020 Shirokov Sergey
#ifndef MODULES_TASK_3_SHIROKOV_S_MULTIPLICATION_MATRIX_MULTIPLICATION_MATRIX_H_
#define MODULES_TASK_3_SHIROKOV_S_MULTIPLICATION_MATRIX_MULTIPLICATION_MATRIX_H_

#include <vector>

#define ERR_INCORRECT_SIZE -2

class SparseMatrix {
public:
    std::vector<double> m_val;
    std::vector<int> m_col_ind;
    std::vector<int> m_row_ptr;
    int m_columns, m_rows, m_non_zero;

    SparseMatrix() : m_columns(0), m_rows(0), m_non_zero(0) {}

    SparseMatrix(const std::vector<double>& t_val,
        const std::vector<int>& t_col_ind,
        const std::vector<int>& t_row_ptr, const int& t_columns,
        const int& t_rows)
        : m_val(t_val),
        m_col_ind(t_col_ind),
        m_row_ptr(t_row_ptr),
        m_columns(t_columns),
        m_rows(t_rows),
        m_non_zero(t_val.size()) {}

    SparseMatrix(const int& t_columns, const int& t_rows)
        : m_columns(t_columns), m_rows(t_rows) {
        m_row_ptr.resize(t_rows + 1);
    }

    SparseMatrix(const std::vector<double>& A, const int& t_columns,
        const int& t_rows);

    const std::vector<double> makeVector() const;
    const std::vector<double> makeColumn(const int& t_col);

    friend const std::vector<double> operator*(const SparseMatrix& A,
        const SparseMatrix& B);
    friend const std::vector<double> operator*(const SparseMatrix& A,
        const std::vector<double>& B);
};


std::vector<double> multiply(SparseMatrix* A, SparseMatrix* B);
std::vector<double> randomSparseMatrix(const int& t_columns, const int& t_rows,
    const int& t_left_bound,
    const int& t_right_bound);

#endif  // MODULES_TASK_3_SHIROKOV_S_MULTIPLICATION_MATRIX_MULTIPLICATION_MATRIX_H_
