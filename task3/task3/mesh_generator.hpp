// ================  MESH_GENERATOR.HPP ================
#pragma once
#ifndef MESH_GENERATOR_HPP
#define MESH_GENERATOR_HPP

#include "mesh.hpp"
#include "area.hpp"

// Абстрактный генератор сетки
class grid_generator
{
protected:
   virtual void read_data(std::string path = directory + "input_data\\") = 0;
   virtual void generate_nodes() = 0;
};

// Генератор двумерной равномерной / неравномерной сетки
class space_grid_generator : public grid_generator
{
public:
   space_grid_generator();

   // Построить сетку и наложить 1-е краевые условия
   void build_mesh(space_grid *&grid, function2D &us, function2D &uc);

private:
   dvector  x, y;       // узлы по каждой оси
   area *area_xy;
   uint32_t nx, ny;
   double   kx, ky;
   mesh_type type;
   uint8_t   nested;    // уровень вложенности (0 — исходная, 1,2,3 — измельчённые)

   void read_data(std::string path = directory + "input_data\\") override;
   void generate_nodes() override;
   void make_bc(space_grid *&grid, function2D &us, function2D &uc);
};

#endif