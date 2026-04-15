// ================  SOLVER.HPP ================
#pragma once
#ifndef SOLVER_HPP
#define SOLVER_HPP

#include "matrix.hpp"
#include "decomposer.hpp"
#include "utilities.hpp"

// Только два метода: прямой (LU) и итерационный (LOS+LU-предобусл.)
enum class method
{
   LU,
   LOS_LU
};

struct result
{
   double residual = 0;
   uint32_t iters = 0;
   std::chrono::microseconds time = std::chrono::microseconds(0);
};

class solver
{
public:
   // Прямой метод (LU-разложение профильной матрицы)
   result solve_by_LU(profile_matrix &A, const dvector &b, dvector &x);

   // Итерационный метод (LOS с LU-предобусловливанием)
   result solve_iterative(
      sparse_matrix &A,
      dvector &b,
      dvector &x,
      method method,
      uint32_t max_iter,
      double eps
   );

private:
   uint32_t max_iter;
   double   eps;

   // Прямой и обратный ход LU-предобусловливателя
   void LU_direct(profile_matrix &A, const dvector &b, dvector &x);
   void LU_reverse(profile_matrix &A, const dvector &b, dvector &x);

   // Локально-оптимальная схема (LOS) + LU
   result LOS_LU(sparse_matrix &A, dvector &b, dvector &x);
};

#endif