#pragma once
#ifndef DIAGONAL_MATRIX_H
#define DIAGONAL_MATRIX_H

#include "head.h"

class diagonal_matrix
{
   friend class fdm;

public:
   diagonal_matrix(const uint32_t _dim, const uint32_t _diags_count, const uint32_t _zero_diags);

   inline uint32_t get_size(void) const { return dim; }

   inline uint32_t get_zero_diags_count() { return zero_diags; }

   double diagonal(const uint32_t di, const uint32_t elem) { return diags[di][elem]; }

   void set_elem(const uint32_t di, const uint32_t elem, const double value) { diags[di][elem] = value; }

   std::unique_ptr<std::vector<double>> dot(const std::vector<double> &vector);

   double dot(const uint32_t row, const std::vector<double> &vector);

   void to_dense();

private:
   uint32_t dim;								// Размер матрицы
   uint32_t zero_diags;						// Кол-во нулевых диагоналей
   std::vector<std::vector<double>> diags;		// Диагонали
   std::vector<int> offset;					// Смещения побочных диагоналей
   // относительно главной
   std::vector<std::vector<double>> matrix;	// Матрица в плотном формате
};

#endif