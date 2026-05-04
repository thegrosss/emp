// ================ FINITE_ELEMENT.HPP ================
#pragma once
#ifndef FINITE_ELEMENT_HPP
#define FINITE_ELEMENT_HPP
#include "head.hpp"

// Конечный элемент ---------------------------------------------
struct finite_elem
{
   std::array<uint32_t, 8> nodes;

   double lambda, sigma, hi;

   uint32_t operator[](const uint32_t index) { return nodes[index]; }

   finite_elem(std::array<uint32_t, 8> _nodes = {}, double _lambda = 0, double _sigma = 0, double _hi = 0)
   {
      nodes = _nodes;
      lambda = _lambda;
      sigma = _sigma;
      hi = _hi;
   }
};

// Грань --------------------------------------------------------
struct face
{
   std::array<uint32_t, 4> nodes;
   uint32_t number;

   face(std::array<uint32_t, 4> _nodes = {}, uint32_t _number = 0)
   {
      nodes = _nodes;
      number = _number;
   }
};

#endif