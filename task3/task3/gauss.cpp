// ================  GAUSS.HPP ================
#include "gauss.hpp"

gauss::gauss()
{
   gauss_points = { 0.0, sqrt(3.0 / 5.0), -sqrt(3.0 / 5.0) };
   gauss_weights = { 8.0 / 9.0, 5.0 / 9.0, 5.0 / 9.0 };
}

double gauss::integrate3D(function3D &f, omega &omega)
{
   double qi, qj, qk;
   double pi, pj, pk;

   double xk = omega.start_point.x;
   double yk = omega.start_point.y;
   double zk = omega.start_point.z;
   double xk1 = omega.end_point.x;
   double yk1 = omega.end_point.y;
   double zk1 = omega.end_point.z;

   double hx = abs(xk1 - xk);
   double hy = abs(yk1 - yk);
   double hz = abs(zk1 - zk);

   double sum = 0.0;

   for (uint8_t i = 0; i < 3; i++)
   {
      qi = gauss_weights[i];
      pi = (xk + xk1 + gauss_points[i] * hx) / 2.0;

      for (uint8_t j = 0; j < 3; j++)
      {
         qj = gauss_weights[j];
         pj = (yk + yk1 + gauss_points[j] * hy) / 2.0;

         for (uint8_t k = 0; k < 3; k++)
         {
            qk = gauss_weights[k];
            pk = (zk + zk1 + gauss_points[k] * hz) / 2.0;

            sum += qi * qj * qk * f(pi, pj, pk);
         }
      }
   }

   return sum * hx * hy * hz / 8.0;
}