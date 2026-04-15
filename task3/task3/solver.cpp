// ================  SOLVER.CPP ================
#include "solver.hpp"

// ===================== Прямой метод (LU-разложение) ========================
result solver::solve_by_LU(profile_matrix &A, const dvector &b, dvector &x)
{
   result res;

   decomposer::LU(A);

   x.resize(b.size());
   double sum = 0.0;

   auto t0 = std::chrono::high_resolution_clock::now();

   // Прямой ход (Ly = b)
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

   // Обратный ход (Ux = y)
   for (int i = static_cast<int>(A.size()) - 1; i >= 0; i--)
   {
      uint32_t i0 = A.ig[i];
      uint32_t i1 = A.ig[i + 1];
      uint32_t j = i - (i1 - i0);

      for (uint32_t k = i0; k < i1; k++)
         x[j++] -= A.ggu[k] * x[i];
   }

   auto t1 = std::chrono::high_resolution_clock::now();

   res.iters = 0;
   res.residual = 0;
   res.time = std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0);
   return res;
}

// ===================== Итерационный диспетчер ===============================
result solver::solve_iterative(sparse_matrix &A, dvector &b, dvector &x,
   method method, uint32_t _max_iter, double _eps)
{
   max_iter = _max_iter;
   eps = _eps;
   return LOS_LU(A, b, x);
}

// ===================== Вспомогательные ходы LU ==============================

void solver::LU_direct(profile_matrix &A, const dvector &b, dvector &x)
{
   x = b;
   for (uint32_t i = 0; i < A.size(); ++i)
   {
      uint32_t i0 = A.ig[i];
      uint32_t i1 = A.ig[i + 1];
      uint32_t j = i - (i1 - i0);
      double sum = 0.0;
      for (uint32_t k = i0; k < i1; ++k)
         sum += A.ggl[k] * x[j++];
      x[i] = (x[i] - sum) / A.di[i];
   }
}

void solver::LU_reverse(profile_matrix &A, const dvector &b, dvector &x)
{
   x = b;
   for (int i = (int)A.size() - 1; i >= 0; --i)
   {
      uint32_t i0 = A.ig[i];
      uint32_t i1 = A.ig[i + 1];
      uint32_t j = i - (i1 - i0);
      for (uint32_t k = i0; k < i1; ++k)
         x[j++] -= A.ggu[k] * x[i];
   }
}

// ===================== LOS + LU-предобусловливание ==========================
/*
  Локально-оптимальная схема (минимизация невязки на каждом шаге):
    r^k, z^k, p^k  — рабочие векторы в пространстве L^{-1}A U^{-1}
*/
result solver::LOS_LU(sparse_matrix &_A, dvector &b, dvector &x)
{
   // Преобразуем в профильный формат – он учитывает fill‑in
   std::unique_ptr<profile_matrix> A_ptr(_A.to_profile());
   profile_matrix &A = *A_ptr;

   // LU-разложение профильной матрицы
   decomposer::LU(A);

   // Вспомогательные векторы (размерность та же)
   uint32_t dim = _A.size();
   dvector r(dim), z(dim), p(dim), LAU(dim), U(dim);
   double alpha, beta;

   x.assign(dim, 0.0);
   auto t0 = std::chrono::high_resolution_clock::now();

   // Прямой и обратный ходы для профильной матрицы (пишем новые функции)
   LU_direct(A, b - *_A.dot(x), r);
   LU_reverse(A, r, z);
   LU_direct(A, *_A.dot(z), p);

   uint32_t k = 0;
   double rr = r * r;
   double eps_sq = eps * eps; // сравнение квадратов для устойчивости

   for (; k < max_iter && rr >= eps_sq; ++k)
   {
      double pp = p * p;
      if (pp == 0.0) break; // защита от вырождения

      alpha = (p * r) / pp;
      x = x + alpha * z;
      r = r - alpha * p;
      rr = r * r;

      // LAU = A^{-T} * r   (обратный ход, затем прямой)
      LU_reverse(A, r, LAU);
      LAU = *_A.dot(LAU);
      LU_direct(A, LAU, LAU);

      beta = -(p * LAU) / pp;
      LU_reverse(A, r, U);
      z = U + beta * z;
      p = LAU + beta * p;
   }

   auto t1 = std::chrono::high_resolution_clock::now();
   return { sqrt(rr), k, std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0) };
}