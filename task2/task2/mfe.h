#pragma once
#ifndef MFE_H
#define MFE_H

#include "mesh.h"
#include "matrix.h"
#include "solver.h"
#include "one_dimensional_search.h"

class mfe
{
   friend class solver;

public:
   enum class method
   {
      SIMPLE_ITERATION,
      NEWTON
   };

   mfe(mesh &mesh);

   void set_functions(const function &_u, const function_f &_f, const function_lambda &_lambda);

   void assembly_global_slae(mesh &mesh, const double t, const double delta_t);

   void initial_condition(mesh &mesh);
   void first_boundary_condition(mesh &mesh, const double t, const double delta_t);

   std::pair<uint32_t, double> solve(mesh &mesh, uint32_t max_iter, double eps);

private:
   std::vector<std::vector<double>> local_A;	// Локальная матрица
   std::vector<double> local_f;			// Локальный вектор

   matrix *global_A;			// Глобальная матрица
   std::vector<double> global_f;	// Глобальный вектор

   matrix *A;				// Глобальная матрица для расчета невязки, если 
                        // используем метод Ньютона
   std::vector<double> nl_f;

   std::vector<double> q;		    // Вектор весов на текущем временном слое
   std::vector<double> q_prev;	    // Вектор весов на предыдущем временном слое
   std::vector<double> q_prev_iter;	// Вектор весов на предыдущем временном слое (для минимизации)

   std::vector<double> exact;	// Точное значение

   function u;				// Неизвестная функция
   function_f f;			// Функция правой части
   function_lambda lambda;		// Функция lambda(du/dx)

   method method_to_solve;
   bool use_relax;

   double dlambda(double q2, double q1, double h, double x, uint32_t var);

   void build_local_matrix(const finite_elem &elem, const double delta_t);
   void build_local_vector(const finite_elem &elem, const double t, const double delta_t);

   void add_relax(const double w);

   double step();
   double residual();
   double error(mesh &mesh);

   void save(std::string dir, std::pair<uint32_t, double> &res);

   void linearization_newton(const finite_elem &elem, const double delta_t);
};
#endif
