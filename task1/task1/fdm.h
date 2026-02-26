#pragma once
#ifndef FDM_H
#define FDM_H

#include "diagonal_matrix.h"
#include "mesh.h"
#include "solver.h"

class fdm
{
public:
   void mesh_to_slae(mesh &mesh, func2D_u &u, func2D_f &f);

   std::pair<uint32_t, double> calculate(std::vector<double> &u);

public:
   diagonal_matrix *A;			// Матрица СЛАУ
   std::vector<double> b;		// Вектор правой части

   std::vector<double> q;

   std::vector<double> exact;	// Точное решение

   double beta = 4;
   double u_beta = 0;

   uint32_t slae_size;			// Размер СЛАУ

   double du_dx(func2D_u &f, point &p, double hx);
   double du_dy(func2D_u &f, point &p, double hy);
   double residual();
};

#endif