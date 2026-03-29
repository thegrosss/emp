#include "mesh.h"

void mesh::save(const std::string dir)
{
   //--------------- КОНЕЧНЫЕ ЭЛЕМЕНТЫ ----------------------------
   std::ofstream fe(dir + "fe.txt");

   if (fe.is_open())
   {
      fe << elems[0].x << std::endl;
      for (uint32_t i = 0; i < elems.size(); i++)
         fe << elems[i].x_next << std::endl;

      fe.close();
   }

   //--------------------- ВРЕМЯ ----------------------------------
   std::ofstream time(dir + "time.txt");

   if (time.is_open())
   {
      for (uint32_t i = 0; i < T.size(); i++)
         time << T[i] << std::endl;

      time.close();
   }
}