// ================ BASIS_FUNCTION.HPP ================
#pragma once
#ifndef BASIS_FUNCTION_HPP
#define BASIS_FUNCTION_HPP

#include "finite_element.hpp"
#include "area.hpp"

/*
  Нумерация узлов внутри элемента:
    2 --- 3        (xk,yk1) --- (xk1,yk1)
    |     |    =>       |             |
    0 --- 1        (xk,yk ) --- (xk1,yk )

  Билинейные базисные функции:
    psi_0 = (xk1-x)(yk1-y) / (hx*hy)
    psi_1 = (x-xk )(yk1-y) / (hx*hy)
    psi_2 = (xk1-x)(y-yk ) / (hx*hy)
    psi_3 = (x-xk )(y-yk ) / (hx*hy)
*/
struct basis_function
{
   double xk, xk1;
   double yk, yk1;
   double hx, hy;

   basis_function(omega &_area)
   {
      xk = _area.start_point.x;  xk1 = _area.end_point.x;
      yk = _area.start_point.y;  yk1 = _area.end_point.y;
      hx = xk1 - xk;
      hy = yk1 - yk;
   }

   // Значение i-й базисной функции в точке point
   double psi(uint32_t ifunc, point2D &point)
   {
      double x = point.x, y = point.y;
      switch (ifunc)
      {
      case 0: return (xk1 - x) * (yk1 - y) / (hx * hy);
      case 1: return (x - xk) * (yk1 - y) / (hx * hy);
      case 2: return (xk1 - x) * (y - yk) / (hx * hy);
      case 3: return (x - xk) * (y - yk) / (hx * hy);
      }
      return 0.0;
   }

   // Частная производная i-й базисной функции:
   //   ivar == 1  =>  d/dx
   //   ivar == 2  =>  d/dy
   double d_psi(uint8_t ifunc, uint8_t ivar, point2D &point)
   {
      double x = point.x, y = point.y;

      if (ivar == 1) // d/dx (аналитически)
      {
         switch (ifunc)
         {
         case 0: return -(yk1 - y) / (hx * hy);
         case 1: return  (yk1 - y) / (hx * hy);
         case 2: return -(y - yk) / (hx * hy);
         case 3: return  (y - yk) / (hx * hy);
         }
      }
      else           // d/dy (аналитически)
      {
         switch (ifunc)
         {
         case 0: return -(xk1 - x) / (hx * hy);
         case 1: return -(x - xk) / (hx * hy);
         case 2: return  (xk1 - x) / (hx * hy);
         case 3: return  (x - xk) / (hx * hy);
         }
      }
      return 0.0;
   }
};

#endif