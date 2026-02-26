#pragma once
#include <vector>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <iomanip>
#include <functional>
#include <string>

typedef std::function<double(double, double)> func2D_u;
typedef std::function<double(double, double, double, double)> func2D_f;

struct point
{
   double x;
   double y;

   point(double _x = 0.0, double _y = 0.0) :
      x(_x), y(_y) { };
};