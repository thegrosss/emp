#pragma once
#ifndef TESTING_MODULE
#define TESTING_MODULE

#include "mesh_generator.h"
#include "mfe.h"

class testing_module
{
public:
   testing_module() { set_functions(); }

   void run();

private:
   void set_functions();

   void save(std::ofstream &file, std::pair<uint32_t, double> &res);

   std::vector<function> u;
   std::vector<function_lambda> lambda;
   std::vector<std::vector<function_f>> f;

   std::vector<std::string> u_names;
   std::vector<std::string> lambda_names;
};

#endif