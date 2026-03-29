#include "solver.h"

// LU разложение матрицы
void solver::LU(matrix &A)
{
   for (uint32_t i = 1; i < A.get_size(); i++)
   {
      A.tape[2][i - 1] /= A.tape[0][i - 1];
      A.tape[0][i] -= A.tape[1][i - 1] * A.tape[2][i - 1];
   }
}

// Решить СЛАУ с помощью LU разложения
void solver::solve_by_LU(matrix A, std::vector<double> &f, std::vector<double> &q)
{
   LU(A);

   // Прямой ход
   q[0] = f[0] / A.tape[0][0];

   // Прямой ход
   for (uint32_t i = 1; i < A.get_size(); i++)
   {
      q[i] = f[i] - q[i - 1] * A.tape[1][i - 1];
      q[i] /= A.tape[0][i];
   }

   // Обратный ход
   for (int32_t i = A.get_size() - 2; i >= 0; i--)
      q[i] -= A.tape[2][i] * q[i + 1];
}