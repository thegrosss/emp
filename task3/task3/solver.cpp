// ================  SOLVER.CPP ================
#include "solver.hpp"

// ======================================= ПРЯМОЙ РЕШАТЕЛЬ ===========================================
result solver::solve_by_LU(profile_matrix &A, const dvector &b, dvector &x)
{
   result res;

   // ----------------------------------------------------------
   // Делаем LU разложение матрицы
   decomposer::LU(A);

   // ----------------------------------------------------------
   // Прямой ход
   x.resize(b.size());

   double sum = 0.0;

   auto time_start = std::chrono::high_resolution_clock::now();
   for (uint32_t i = 0; i < A.size(); i++)
   {
      uint32_t i0 = A.ig[i];
      uint32_t i1 = A.ig[i + 1];

      uint32_t j = i - (i1 - i0);

      sum = 0.0;

      for (uint32_t k = i0; k < i1; k++)
         sum += A.ggl[k] * x[j++];

      x[i] = (b[i] - sum) / A.di[i];
   }

   // ----------------------------------------------------------
   // Обратный ход
   for (int i = A.size() - 1; i >= 0; i--)
   {
      uint32_t i0 = A.ig[i];
      uint32_t i1 = A.ig[i + 1];

      uint32_t j = i - (i1 - i0);

      for (uint32_t k = i0; k < i1; k++)
         x[j++] -= A.ggu[k] * x[i];
   }

   auto time_end = std::chrono::high_resolution_clock::now();

   res.iters = 0;
   res.residual = 0;
   res.time = std::chrono::duration_cast<std::chrono::microseconds>(time_end - time_start);

   return res;
}


// =================================== ИТЕРАЦИОННЫЕ РЕШАТЕЛИ =========================================
result solver::solve_iterative(sparse_matrix &A, dvector &b, dvector &x, method method, uint32_t _max_iter, double _eps)
{
   max_iter = _max_iter;
   eps = _eps;

   if (method == method::LOS_LU)
      return LOS_LU(A, b, x);

   return BCGSTAB_LU(A, b, x);
}


void solver::LU_direct(sparse_matrix &A, const dvector &b, dvector &x)
{
   x = b;

   for (uint32_t i = 0; i < x.size(); i++)
   {
      double sum = 0.0;

      for (uint32_t j = A.ig[i]; j < A.ig[i + 1]; j++)
         sum += A.ggl[j] * x[A.jg[j]];

      x[i] -= sum;
      x[i] /= A.di[i];
   }
}

void solver::LU_reverse(sparse_matrix &A, const dvector &b, dvector &x)
{
   x = b;

   for (__int64 i = A.size() - 1; i >= 0; i--)
   {
      for (uint32_t j = A.ig[i]; j < A.ig[i + 1]; j++)
         x[A.jg[j]] -= A.ggu[j] * x[i];
   }
}

result solver::LOS_LU(sparse_matrix &_A, dvector &b, dvector &x)
{
   sparse_matrix A = _A;

   result res;

   uint32_t dim = _A.size();

   dvector r(dim), z(dim), p(dim);
   dvector LAU(dim), U(dim);

   double alpha, beta;


   x.resize(dim, 0);

   decomposer::LU(A);

   auto time_start = std::chrono::high_resolution_clock::now();

   LU_direct(A, b - *_A.dot(x), r);
   LU_reverse(A, r, z);
   LU_direct(A, *_A.dot(z), p);


   uint32_t k = 0;
   double rr = (r * r);

   for (; k < max_iter && rr >= eps; k++)
   {
      double pp = p * p;
      alpha = (p * r) / pp;

      x = x + alpha * z;

      rr = (r * r);

      r = r - alpha * p;

      LU_reverse(A, r, LAU);
      LAU = *_A.dot(LAU);
      LU_direct(A, LAU, LAU);

      beta = -(p * LAU) / pp;

      LU_reverse(A, r, U);

      z = U + beta * z;

      p = LAU + beta * p;
   }

   auto time_end = std::chrono::high_resolution_clock::now();

   res.iters = k;
   res.residual = rr;
   res.time = std::chrono::duration_cast<std::chrono::microseconds>(time_end - time_start);

   r.clear();
   z.clear();
   p.clear();
   LAU.clear();
   U.clear();

   return res;
}

result solver::BCGSTAB_LU(sparse_matrix &_A, dvector &b, dvector &x)
{
   sparse_matrix A = _A;

   uint32_t dim = A.size();

   result res;

   dvector r0(dim), r(dim), z(dim), p(dim), y(dim);
   dvector LAUz(dim), LAUp(dim);

   double alpha, beta, gamma;

   x.resize(b.size(), 0);

   decomposer::LU(A);

   auto time_start = std::chrono::high_resolution_clock::now();

   LU_direct(A, b - *_A.dot(x), r0);

   z = r0;
   r = r0;


   uint32_t k = 0;
   double rr = (r * r);

   for (; k < max_iter && rr >= eps; k++)
   {
      LU_reverse(A, z, LAUz);
      LAUz = *_A.dot(LAUz);
      LU_direct(A, LAUz, LAUz);

      double r_r0 = (r * r0);

      alpha = r_r0 / (LAUz * r0);

      p = r - alpha * LAUz;

      LU_reverse(A, p, LAUp);
      LAUp = *_A.dot(LAUp);
      LU_direct(A, LAUp, LAUp);

      gamma = (p * LAUp) / (LAUp * LAUp);
      y = y + alpha * z + gamma * p;

      r = p - gamma * LAUp;

      beta = alpha * (r * r0) / (gamma * r_r0);

      z = r + beta * z - beta * gamma * LAUz;

      LU_reverse(A, y, x);

      rr = (r * r);
   }

   auto time_end = std::chrono::high_resolution_clock::now();

   res.iters = k;
   res.residual = rr;
   res.time = std::chrono::duration_cast<std::chrono::microseconds>(time_end - time_start);

   r.clear();
   r0.clear();
   z.clear();
   p.clear();
   y.clear();
   LAUp.clear();
   LAUz.clear();

   return res;
}