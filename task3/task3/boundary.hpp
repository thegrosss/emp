// ================ BOUNDARY.HPP ================
#pragma once
#ifndef BOUNDARY_HPP
#define BOUNDARY_HPP

#include "head.hpp"

enum class boundary_type
{
   DIRICHLET = 1,
   NEUMANN = 2
};

struct boundary_cond
{
   uint8_t edge;
   boundary_type type;
};

struct dirichlet_cond
{
   uint32_t node;
   double value;
};

#endif