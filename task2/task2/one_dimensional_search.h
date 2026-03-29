#pragma once
#ifndef ONE_DIMENSIONAL_SEARCH
#define ONE_DIMENSIONAL_SEARCH

#include "head.h"
#include "matrix.h"

struct interval
{
   double a;
   double b;

   interval(double _a, double _b) :
      a(_a), b(_b) { };
};

class one_dimensional_search
{
public:
   static double minimize(matrix &A, std::vector<double> &q, std::vector<double> &q_prev, std::vector<double> &b, double eps);

   static interval find_interval(min_func &f, double x0, double eps);

   static double golden_ratio(min_func &f, interval &interval, double eps);
};

#endif