// ================  GAUSS.CPP ================
#include "gauss.hpp"

gauss::gauss()
{
   gauss_points = { 0.0, sqrt(3.0 / 5.0), -sqrt(3.0 / 5.0) };
   gauss_weights = { 8.0 / 9.0, 5.0 / 9.0,  5.0 / 9.0 };
}

// Формула: I ≈ (hx*hy/4) * Σ_i Σ_j w_i * w_j * f(xi, yj)
// где xi = (xk + xk1 + t_i * hx) / 2,  yj = (yk + yk1 + t_j * hy) / 2
double gauss::integrate2D(function2D &f, omega &omega)
{
   double xk = omega.start_point.x, xk1 = omega.end_point.x;
   double yk = omega.start_point.y, yk1 = omega.end_point.y;
   double hx = std::abs(xk1 - xk);
   double hy = std::abs(yk1 - yk);

   double sum = 0.0;

   for (uint8_t i = 0; i < 3; i++)
   {
      double qi = gauss_weights[i];
      double pi = (xk + xk1 + gauss_points[i] * hx) / 2.0;

      for (uint8_t j = 0; j < 3; j++)
      {
         double qj = gauss_weights[j];
         double pj = (yk + yk1 + gauss_points[j] * hy) / 2.0;

         sum += qi * qj * f(pi, pj);
      }
   }

   return sum * hx * hy / 4.0;
}