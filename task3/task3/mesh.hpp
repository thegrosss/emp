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

// ===================== Двумерная расчётная сетка ==================
class space_grid : public grid
{
   friend class space_grid_generator;

public:
   // Сетка nx * ny элементов, (nx+1)*(ny+1) узлов
   space_grid(const uint32_t _nx = 0, const uint32_t _ny = 0)
   {
      nx = _nx;
      ny = _ny;
      elems.resize(nx * ny);
      points.resize((nx + 1) * (ny + 1));
      edges.resize(2 * (nx + ny)); // граничные рёбра
   }

   ~space_grid()
   {
      elems.clear();
      points.clear();
      edges.clear();
      dirichlet.clear();
   }

   inline uint32_t get_width(void)       const { return nx; }
   inline uint32_t get_height(void)      const { return ny; }
   inline uint32_t get_nodes_count(void) const override { return points.size(); }
   inline uint32_t get_elems_count(void) const { return elems.size(); }
   inline uint32_t get_edges_count(void) const { return edges.size(); }

   inline finite_elem &get_elem(const uint32_t index) { return elems[index]; }
   inline point2D &get_point(const uint32_t index) { return points[index]; }

   inline uint32_t get_first_bound_count(void) const { return dirichlet.size(); }
   inline dirichlet_cond &get_dirichlet_cond(const uint32_t num) { return dirichlet[num]; }

   void save(std::string path = directory + "mesh\\") override;

   // Установка параметров по всем элементам
   void set_lambda(const double lambda)
   {
      for (auto &e : elems) e.lambda = lambda;
   }
   void set_sigma(const double sigma)
   {
      for (auto &e : elems) e.sigma = sigma;
   }
   void set_hi(const double hi)
   {
      for (auto &e : elems) e.hi = hi;
   }

private:
   uint32_t nx, ny;

   std::vector<point2D>       points;
   std::vector<finite_elem>   elems;
   std::vector<edge>          edges;    // граничные рёбра (1-е краевое условие)

   std::vector<dirichlet_cond> dirichlet;
};

#endif