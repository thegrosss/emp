// ================ UTILITIES.HPP ================
#pragma once
#ifndef UTILITIES_HPP
#define UTILITIES_HPP

#include "head.hpp"

inline static double norm(const dvector &vector)
{
   double scalar = 0.0;

   for (const auto &it : vector)
      scalar += it * it;

   return sqrt(scalar);
}

inline static dvector operator+(const dvector &v1, const dvector &v2)
{
   dvector result = v1;

   for (uint32_t i = 0; i < result.size(); i++)
      result[i] += v2[i];

   return result;
}

inline static dvector operator-(const dvector &v1, const dvector &v2)
{
   dvector result = v1;

   for (uint32_t i = 0; i < result.size(); i++)
      result[i] -= v2[i];

   return result;
}

inline static double operator*(const dvector &v1, const dvector &v2)
{
   double scalar = 0.0;

   for (uint32_t i = 0; i < v1.size(); i++)
      scalar += v1[i] * v2[i];

   return scalar;
}

inline static dvector operator*(const double c, const dvector &v)
{
   dvector result = v;

   for (uint32_t i = 0; i < result.size(); i++)
      result[i] *= c;

   return result;
}

#endif