// ================  DECOMPOSER.HPP ================
#include "decomposer.hpp"

void decomposer::LU(sparse_matrix &A)
{
   for (uint32_t i = 0; i < A.size(); i++)
   {
      double sum_di = 0.0;

      uint32_t i0 = A.ig[i];
      uint32_t i1 = A.ig[i + 1];

      for (uint32_t k = i0; k < i1; k++)
      {
         uint32_t j = A.jg[k];

         uint32_t j0 = A.ig[j];
         uint32_t j1 = A.ig[j + 1];

         uint32_t ik = i0;
         uint32_t kj = j0;

         double sum_l = 0.0;
         double sum_u = 0.0;

         while (ik < k && kj < j1)
         {
            if (A.jg[ik] == A.jg[kj])
            {
               sum_l += A.ggl[ik] * A.ggu[kj];
               sum_u += A.ggu[ik] * A.ggl[kj];
               ik++;
               kj++;
            }
            else
               A.jg[ik] > A.jg[kj] ? kj++ : ik++;
         }

         A.ggl[k] -= sum_l;
         A.ggu[k] = (A.ggu[k] - sum_u) / A.di[j];
         sum_di += A.ggl[k] * A.ggu[k];
      }

      A.di[i] -= sum_di;
   }
}

void decomposer::LU(profile_matrix &A)
{
   for (uint32_t i = 1; i < A.size(); i++)
   {
      double sum_di = 0.0;

      uint32_t j0 = i - (A.ig[i + 1] - A.ig[i]);

      for (uint32_t ii = A.ig[i]; ii < A.ig[i + 1]; ii++)
      {
         uint32_t j = ii - A.ig[i] + j0;
         uint32_t jbeg = A.ig[j];
         uint32_t jend = A.ig[j + 1];

         if (jbeg < jend)
         {
            uint32_t i0 = j - (jend - jbeg);
            uint32_t jjbeg = std::max(j0, i0);
            uint32_t jjend = std::min(j, i - 1);

            double sum_l = 0.0;
            double sum_u = 0.0;

            for (uint32_t k = 0; k < jjend - jjbeg; k++)
            {
               uint32_t indau = A.ig[j] + jjbeg - i0 + k;
               uint32_t indal = A.ig[i] + jjbeg - j0 + k;
               sum_l += A.ggu[indau] * A.ggl[indal];
            }
            A.ggl[ii] -= sum_l;

            for (uint32_t k = 0; k < jjend - jjbeg; k++)
            {
               uint32_t indal = A.ig[j] + jjbeg - i0 + k;
               uint32_t indau = A.ig[i] + jjbeg - j0 + k;
               sum_u += A.ggu[indau] * A.ggl[indal];
            }
            A.ggu[ii] -= sum_u;
         }
         A.ggu[ii] /= A.di[j];

         sum_di += A.ggl[ii] * A.ggu[ii];
      }

      A.di[i] -= sum_di;
   }
}