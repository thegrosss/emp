// ================  GAUSS.HPP ================
#pragma once
#ifndef GAUSS_HPP
#define GAUSS_HPP

#include "head.hpp"
#include "area.hpp"

class gauss
{
public:
   gauss();

   // Двумерное интегрирование квадратурами Гаусса (3x3 точки)
   double integrate2D(function2D &f, omega &omega);

private:
   std::array<double, 3> gauss_points;
   std::array<double, 3> gauss_weights;
};

#endif