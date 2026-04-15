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
   ~mfe() { delete A; }

   void assembly_global_matrix_and_vector(r_function2D &fs, r_function2D &fc);

   void add_dirichlet();

   inline void set_w(const double _w) { w = _w; }

   result solve(method method_to_solve);

public:
   double w;

   // Локальные матрицы 4×4 (билинейные базисные функции)
   std::vector<dvector> local_p; // λ*G - ω²χ*M  (блок P_ij)
   std::vector<dvector> local_c; // ω*σ*M         (блок C_ij)
   std::vector<dvector> G;       // матрица жёсткости
   std::vector<dvector> M;       // матрица масс

   dvector local_vec_fs; // правая часть (sin-часть)
   dvector local_vec_fc; // правая часть (cos-часть)

private:
   void build_local_matrices(finite_elem &elem);
   void build_local_vectors(finite_elem &elem, r_function2D &fs, r_function2D &fc);

   void add_to_global_matrix(const uint32_t i, const uint32_t j, const double val);

   space_grid *grid;  // расчётная сетка
   sparse_matrix *A;     // глобальная матрица
   gauss          gauss; // интегрирование по Гауссу

public:
   dvector q; // вектор решения
   dvector b; // вектор правой части
};

#endif