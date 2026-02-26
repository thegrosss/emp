#include "solver.h"

// Задать параметры для решения СЛАУ
void solver::set_parameters(const uint32_t _max_iter,
   const double _eps)
{
   max_iter = _max_iter;
   eps = _eps;
}

// Задать начальное приближение
void solver::set_initial_approx(std::vector<double> &vector, const uint32_t size)
{
   vector.clear();
   vector.resize(size, 0);
}

// Решить систему
std::pair<uint32_t, double> solver::solve(double w,
   diagonal_matrix &A,
   std::vector<double> &b,
   std::vector<double> &result)
{
   set_initial_approx(result, A.get_size());

   norm_b = norm(b);

   uint32_t iters = 0;

   do
   {
      iters++;
      Gauss_Seidel(w, A, b, result);

   } while (iters < max_iter && residual(A, b, result) >= eps);

   return std::make_pair(iters, residual(A, b, result));
}

// Метод Гаусса - Зейделя
void solver::Gauss_Seidel(double w,
   diagonal_matrix &A,
   std::vector<double> &b,
   std::vector<double> &result)
{
   for (uint32_t i = 0; i < result.size(); i++)
      result[i] += w * (b[i] - A.dot(i, result)) / A.diagonal(0, i);
}

// Посчитать норму вектора
double solver::norm(const std::vector<double> &vector)
{
   double scalar = 0.0;

   for (const auto &it : vector)
      scalar += it * it;

   return sqrt(scalar);
}

// Посчитать относительную невязку
double solver::residual(diagonal_matrix &A,
   const std::vector<double> &b,
   std::vector<double> &result)
{
   auto tmp_vector = A.dot(result);
   for (uint32_t i = 0; i < result.size(); i++)
      (*tmp_vector)[i] -= b[i];

   return norm(*tmp_vector) / norm_b;
}