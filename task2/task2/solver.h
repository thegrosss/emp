#pragma once
#ifndef SOLVER_H
#define SOLVER_H

#include "head.h"
#include "matrix.h"

class solver
{
public:
   void LU(matrix &A);

   void solve_by_LU(matrix A, std::vector<double> &f, std::vector<double> &q);


private:
   uint32_t max_iter;
   double eps;
};

#endif