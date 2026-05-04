// ================  MFE.HPP ================
#pragma once
#ifndef MFE_HPP
#define MFE_HPP

#include "mesh.hpp"
#include "matrix.hpp"
#include "basis_function.hpp"
#include "gauss.hpp"
#include "portrait_generator.hpp"
#include "solver.hpp"
#include "utilities.hpp"

class mfe
{
public:
   mfe(space_grid &mesh);

   void assembly_global_matrix_and_vector(r_function3D &fs, r_function3D &fc);

   void add_dirichlet();

   inline void set_w(const double _w) { w = _w; }

   result solve(method method_to_solve);

public:
   double w;

   std::vector<dvector> local_p;
   std::vector<dvector> local_c;

   std::vector<dvector> G;
   std::vector<dvector> M;

   dvector local_vec_fs;
   dvector local_vec_fc;

   void build_local_matrices(finite_elem &elem);
   void build_local_vectors(finite_elem &elem, r_function3D &fs, r_function3D &fc);

   void add_to_global_matrix(const uint32_t i, const uint32_t j, const double val);

   space_grid *grid;		// Сетка

   sparse_matrix *A;		// Глобальная матрица
   dvector q;	// Вектор весов
   dvector b;	// Вектор правой части

   gauss gauss;			// Гаусс (для интегрирвоания)
};

#endif