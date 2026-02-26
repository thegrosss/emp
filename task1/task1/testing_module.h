#pragma once
#ifndef TESTING_MODULE_H
#define TESTING_MODULE_H

#include "head.h"
#include "mesh_generator.h"
#include "fdm.h"

class testing_module
{
public:
   void set_functions();

   void run_tests();

private:
   std::vector<std::pair<func2D_u, func2D_f>> functions;

   std::vector<double> exact;

   void calc_exact(const func2D_u &u, mesh &mesh);

   double error(const std::vector<double> &u);

   void save(uint32_t test_num, std::pair<uint32_t, double> &result, const std::vector<double> &u);
};

#endif
