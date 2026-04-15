// ================  MATRIX.HPP ================
#pragma once
#ifndef MATRIX_HPP
#define MATRIX_HPP

#include "head.hpp"

enum class matrix_type
{
   DENSE,
   PROFILE,
   SPARSE
};

class profile_matrix
{
   friend class solver;
   friend class decomposer;
   friend class sparse_matrix;

public:
   profile_matrix(uint32_t _size) { dim = _size; }

   std::unique_ptr<dvector> dot(const dvector &vector);

   inline uint32_t size() const { return dim; };

   void to_dense();

   void save(matrix_type type, std::string path = directory + "matrix\\");

public:
   uint32_t dim;

   std::vector<uint32_t> ig;
   dvector di, ggl, ggu;

   std::vector<dvector> dense_matrix;
};


class sparse_matrix
{
   friend class solver;
   friend class mfe;
   friend class decomposer;

public:
   sparse_matrix(uint32_t _size) { dim = _size; }

   std::unique_ptr<dvector> dot(const dvector &vector);

   inline uint32_t size() const { return dim; };

   void to_dense();
   profile_matrix *to_profile();

   void save(matrix_type type, std::string path = directory + "matrix\\");

public:
   uint32_t dim;

   std::vector<uint32_t> ig, jg;
   dvector di, ggl, ggu;

   std::vector<dvector> dense_matrix;
};

#endif