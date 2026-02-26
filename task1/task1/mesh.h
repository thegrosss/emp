#pragma once
#ifndef MESH_H
#define MESH_H

#include "area.h"

// Узел сетки: 
// координаты, 
// позиции в массивах X и Y, 
// тип (внутренний, граничный, фиктивный)
// тип краевого
// коэффициенты лямбда и гамма
struct node
{
   point p;
   uint32_t x_pos, y_pos;
   int type;				// -1 - фиктивный
   // 0 - внутренний
   // 1,2,3,4 - граничный
   border::bound_cond bc;
   double lambda;
   double gamma;

   node(
      point &_p, uint32_t i, uint32_t j,
      int _type = -1,
      border::bound_cond _bc = border::bound_cond::NONE,
      double _lambda = 0.0, double _gamma = 0.0) :
      p(_p), x_pos(i), y_pos(j), type(_type), bc(_bc), lambda(_lambda), gamma(_gamma)
   {
   };
};

class mesh
{
public:
   enum class mesh_type
   {
      UNIFORM,
      NONUNIFORM
   };

   mesh() : width(0), height(0), type(mesh_type::UNIFORM) {}

   inline void add_node(const node &node) { nodes.push_back(node); }

   inline void set_width(const uint32_t _width) { width = _width; }
   inline void set_height(const uint32_t _height) { height = _height; }
   inline void set_type(const mesh_type _type) { type = _type; }

   inline uint32_t get_width(void) const { return width; }
   inline uint32_t get_height(void) const { return height; }
   inline mesh_type get_type(void) const { return type; }

   inline uint32_t size(void) const { return nodes.size(); }

   node &operator[] (const uint32_t index) { return nodes[index]; }

   void save();

   inline void clear() { nodes.clear(); }

private:
   uint32_t width;				// Ширина сетки
   uint32_t height;			// Высота сетки

   mesh_type type;				// Тип сетки

   std::vector<node> nodes;	// Узлы
};

#endif