// ================  TESTS.HPP ================
#pragma once
#ifndef TESTS_HPP
#define TESTS_HPP

#include "mesh.hpp"
#include "mesh_generator.hpp"
#include "mfe.hpp"
#include "solver.hpp"
#include "utilities.hpp"

struct test_functions
{
   function3D us, uc;
   r_function3D fs, fc;
};

class tests
{
public:
   tests();

   void run();

private:
   void omega_tests();

   void lambda_tests();

   void sigma_tests();

   void hi_tests();

   void calc_exact(space_grid &grid);

   dvector exact;

   std::vector<test_functions> funcs;
};

#endif