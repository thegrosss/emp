// ================  DECOMPOSER.HPP ================
#pragma once
#ifndef DECOMPOSER_HPP
#define DECOMPOSER_HPP

#include "matrix.hpp"

class decomposer
{
public:
   static void LU(sparse_matrix &A);

   static void LU(profile_matrix &A);
};

#endif