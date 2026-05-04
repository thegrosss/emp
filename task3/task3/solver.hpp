// ================  SOLVER.HPP ================
#pragma once
#ifndef SOLVER_HPP
#define SOLVER_HPP

#include "matrix.hpp"
#include "decomposer.hpp"
#include "utilities.hpp"

enum class method
{
   LU,
   LOS_LU,
   BCGSTAB_LU
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
   result solve_iterative(
      sparse_matrix &A,
      dvector &b,
      dvector &x,
      method method,
      uint32_t max_iter, double eps
   );

   result solve_by_LU(profile_matrix &A, const dvector &b, dvector &x);

private:
   // Параметры для итерационного решателя
   uint32_t max_iter;
   double eps;

   // Итерационные методы -------------------------------------------------------------
   // Локально - оптимальная схема + LU предобуславливание
   void LU_direct(sparse_matrix &A, const dvector &b, dvector &x);
   void LU_reverse(sparse_matrix &A, const dvector &b, dvector &x);

   result LOS_LU(sparse_matrix &A, dvector &b, dvector &x);

   // Метод бисопряженных градиаентов стабилизированный + LU предобуславливание
   result BCGSTAB_LU(sparse_matrix &A, dvector &b, dvector &x);
};

#endif