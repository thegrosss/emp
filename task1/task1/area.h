#pragma once
#include "head.h"

// Граница: задает пределы (2 точки) и тип краевого
struct border
{
   enum class bound_cond
   {
      NONE,
      DIRICHLET,
      NEUMANN
   };

   point limits[2];
   bound_cond bc = bound_cond::NONE;
};

// Элементарная область: задает 4 границы и коэффициенты
struct area
{
   border borders[4];
   double lambda = 0.0;
   double gamma = 0.0;
};