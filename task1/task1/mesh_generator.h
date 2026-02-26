#pragma once
#ifndef MESH_GENERATOR_H
#define MESH_GENERATOR_H

#include "head.h"
#include "mesh.h"
#include "area.h"

class mesh_generator
{
public:
   mesh_generator()
   {
      type = mesh::mesh_type::UNIFORM;
      n_omega = 0;
      nx = ny = 0;
      nesting = 0;
   }

   void build_mesh(mesh &mesh);

private:
   mesh::mesh_type type;

   uint32_t n_omega;						 // Кол-во подобластей
   uint32_t nx, ny;						 // Кол-во элементов в массивах
                                     // X_lines и Y_lines соответственно

   std::vector<area> areas;				 // Массив с подобластями

   std::vector<double> X_lines;			 // Координатные линии по X
   std::vector<double> Y_lines;			 // Координатные линии по Y

   std::vector<std::vector<uint32_t>> part; // Массив с информацией о разбиениях
   std::vector<std::vector<double>> kr;	  // Массив с коэффициентами разрядки

   uint32_t nesting;          // Вложенность сетки: 0 - обычная
                              //                    1 - вложенная
                              //                    2 - дважды вложенная
                              //                    3 - трижды вложенная

   std::vector<double> X;	   // Массив X с учетом разбиений
   std::vector<double> Y;	   // Массив Y с учетом разбиений

   void input();

   void generate_xy(mesh::mesh_type &type);

   uint32_t what_type(const point &p, const uint32_t i, const uint32_t j);

   uint32_t what_border(const point &p, const uint32_t area_num);

   bool is_in_area(const point &p, uint32_t &area_num);
};

#endif