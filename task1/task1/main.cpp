#include <iostream>
#include "testing_module.h"

uint32_t main()
{
   testing_module test;

   test.set_functions();
   test.run_tests();

   return 0;
}