#include "one_dimensional_search.h"

// Поиск интервала, содержащего минимум функции
interval one_dimensional_search::find_interval(min_func &f, double x0, double eps)
{
   double delta = eps / 2.0;

   double xk_1, xk1;
   double h;

   double f1 = f(x0);
   double f2 = f(x0 + delta);

   if (f1 > f2)
   {
      xk1 = x0 + delta;
      h = delta;
   }
   else
   {
      xk1 = x0 - delta;
      h = -delta;
   }

   f2 = f(xk1);
   do
   {
      xk_1 = x0;
      x0 = xk1;
      f1 = f2;

      h *= 2;
      xk1 = x0 + h;
      f2 = f(xk1);

   } while (f1 > f2);

   double a = xk_1;
   double b = xk1;

   if (b < a)
      std::swap(a, b);

   return interval(a, b);
}

// Поиск минимума методом золотого сечения
double one_dimensional_search::golden_ratio(min_func &f, interval &interval, double eps)
{
   double delta = eps / 2.0;

   double a = interval.a;
   double b = interval.b;

   double x1 = a + (3.0 - sqrt(5.0)) / 2.0 * (b - a);
   double x2 = a + (sqrt(5.0) - 1) / 2.0 * (b - a);

   double f1 = f(x1);
   double f2 = f(x2);

   while (abs(b - a) >= eps)
   {
      if (f1 > f2)
      {
         a = x1;
         x1 = x2;
         f1 = f2;
         x2 = a + (sqrt(5.0) - 1) / 2.0 * (b - a);
         f2 = f(x2);
      }
      else
      {
         b = x2;
         x2 = x1;
         f2 = f1;
         x1 = a + (3.0 - sqrt(5.0)) / 2.0 * (b - a);
         f1 = f(x1);
      }
   }

   return (a + b) / 2.0;
}

// Минимизировать функционал
double one_dimensional_search::minimize(matrix &A, std::vector<double> &q, std::vector<double> &q_prev, std::vector<double> &b, double eps)
{
   std::vector<double> tmp;

   min_func f = [&](double w)
      {
         tmp = q;

         for (uint32_t i = 0; i < q.size(); i++)
            tmp[i] = tmp[i] * w + (1 - w) * q_prev[i];

         auto prod = A.dot(tmp);

         for (uint32_t i = 0; i < q.size(); i++)
            (*prod)[i] -= b[i];

         return norm(*prod);
      };

   interval interval = find_interval(f, -1.0, eps);

   return golden_ratio(f, interval, eps);
}