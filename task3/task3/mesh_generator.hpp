// ================  MESH_GENERATOR.HPP ================
#pragma once
#ifndef MESH_GENERATOR_HPP
#define MESH_GENERATOR_HPP

#include "mesh.hpp"
#include "area.hpp"

// ---------------------------------------------------------------------------
// Абстрактный класс генератора сетки
class grid_generator
{
protected:
   virtual void read_data(std::string path = directory + "input_data\\") = 0;

   virtual void generate_nodes() = 0;
};


// =================== ГЕНЕРАЦИЯ СЕТКИ ПО ПРОСТРАНСТВУ ====================
class space_grid_generator : public grid_generator
{
public:
   space_grid_generator();

   void build_mesh(space_grid *&grid, function3D &us, function3D &uc);

private:
   // -------------------------------------------------------------
   dvector x, y, z;
   area *area_xyz;
   uint32_t nx, ny, nz;
   double kx, ky, kz;
   mesh_type type;
   uint8_t nested;

   std::vector<boundary_cond> bconds;

   // -------------------------------------------------------------
   void read_data(std::string path = directory + "input_data\\") override;

   void generate_nodes() override;

   void make_bc(space_grid *&grid, function3D &us, function3D &uc);
};

#endif