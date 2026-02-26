#pragma once
#ifndef SOLVER_H
#define SOLVER_H

#include "diagonal_matrix.h"

class solver
{
public:
   solver(uint32_t _max_iter, double _eps) :
      max_iter(_max_iter), eps(_eps), norm_b(0.0) { }

   void set_parameters(const uint32_t _max_iter,
      const double _eps);

   std::pair<uint32_t, double> solve(double w,
      diagonal_matrix &A,
      std::vector<double> &b,
      std::vector<double> &result);

private:
   uint32_t max_iter;		// Максимальное число итераций

   double eps;				// Точность решения

   double norm_b;			// Норма вектора правой части(считаем один
   // раз в начале, чтобы не пересчитывать каждый
   // раз, когда нужно будет посчитать невящку)

   void Gauss_Seidel(double w,
      diagonal_matrix &A,
      std::vector<double> &b,
      std::vector<double> &result);

   void set_initial_approx(std::vector<double> &vector, const uint32_t size);

   double residual(diagonal_matrix &A, const std::vector<double> &b, std::vector<double> &result);

   double norm(const std::vector<double> &vector);
};

#endif