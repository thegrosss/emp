// ================  PORTRAIT_GENERATOR.HPP ================
#pragma once
#ifndef PORTRAIT_GENERATOR_HPP
#define PORTRAIT_GENERATOR_HPP

#include "mesh.hpp"

class portrait_generator
{
public: static void portrait(space_grid &grid, std::vector<uint32_t> &ig, std::vector<uint32_t> &jg);
};

#endif