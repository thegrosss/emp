#pragma once
#ifndef MATRIX_H
#define MATRIX_H

#include "head.h"

class matrix
{
   friend class solver;

public:
   matrix(const uint32_t _size, const uint32_t _tape_width);

   inline void set_elem(const uint32_t diag, const uint32_t elem, const double value) { tape[diag][elem] = value; }
   inline void add_to_elem(const uint32_t diag, const uint32_t elem, const double value) { tape[diag][elem] += value; }

   void set_elem_to_null();

   inline uint32_t get_size(void) const { return size; }

   inline uint32_t get_tape_width(void) const { return tape_width; }

   std::unique_ptr<std::vector<double>> dot(const std::vector<double> &vector);

   void to_dense(const std::string dir = directory);

public:
   uint32_t size;
   uint32_t tape_width;

   std::vector<std::vector<double>> tape;

   std::vector<std::vector<double>> mat;
};

#endif
