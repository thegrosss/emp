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

   double integrate3D(function3D &f, omega &omega);

private:
   std::array<double, 3> gauss_points;
   std::array<double, 3> gauss_weights;
};

#endif