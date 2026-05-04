// ================  MESH.HPP ================
#pragma once
#ifndef MESH_HPP
#define MESH_HPP

#include "finite_element.hpp"
#include "boundary.hpp"
#include "area.hpp"

enum class mesh_type
{
   UNIFORM,
   NONUNIFORM
};


// ==================== Абстрактный класс сетки =====================
class grid
{
public:
   grid() : type(mesh_type::UNIFORM) {};

   virtual void set_type(const mesh_type _type) { type = _type; }

   virtual inline uint32_t get_nodes_count(void) const = 0;
   virtual inline mesh_type get_type(void) const { return type; }

   virtual void save(std::string path = directory + "mesh\\") = 0;

protected:
   mesh_type type;
};


// ===================== Сетка по пространству ======================
class space_grid : public grid
{
   friend class space_grid_generator;

public:
   space_grid(const uint32_t _nx = 0, const uint32_t _ny = 0, const uint32_t _nz = 0)
   {
      nx = _nx;
      ny = _ny;
      nz = _nz;

      elems.resize(nx * ny * nz);
      points.resize((nx + 1) * (ny + 1) * (nz + 1));
      faces.resize(2 * (points.size() - elems.size() - nx - ny - nz - 1));
   };

   ~space_grid()
   {
      elems.clear();
      points.clear();
      faces.clear();
      dirichlet.clear();
      //neumann.clear();
   }

   inline uint32_t get_width(void) const { return nx; }
   inline uint32_t get_lendth(void) const { return ny; }
   inline uint32_t get_height(void) const { return nz; }
   inline uint32_t get_nodes_count(void) const override { return points.size(); }
   inline uint32_t get_elems_count(void) const { return elems.size(); }
   inline uint32_t get_faces_count(void) const { return faces.size(); }

   inline finite_elem &get_elem(const uint32_t index) { return elems[index]; }
   inline point3D &get_point(const uint32_t index) { return points[index]; }

   inline uint32_t get_first_bound_count(void) const { return dirichlet.size(); }
   inline dirichlet_cond &get_dirichlet_cond(const uint32_t num) { return dirichlet[num]; }

   void save(std::string path = directory + "mesh\\") override;


   // ДЛЯ ТЕСТОВ
   void set_lambda(const double lambda)
   {
      for (uint32_t i = 0; i < elems.size(); i++)
         elems[i].lambda = lambda;
   }

   void set_sigma(const double sigma)
   {
      for (uint32_t i = 0; i < elems.size(); i++)
         elems[i].sigma = sigma;
   }

   void set_hi(const double hi)
   {
      for (uint32_t i = 0; i < elems.size(); i++)
         elems[i].hi = hi;
   }

private:
   uint32_t nx, ny, nz;

   std::vector<point3D> points;
   std::vector<finite_elem> elems;
   std::vector<face> faces;

   std::vector<dirichlet_cond> dirichlet;
   //std::vector<neumann_cond> neumann;
};

#endif