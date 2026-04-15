// ================ AREA.HPP ================
#pragma once
#ifndef AREA_HPP
#define AREA_HPP

// Точка в двумерном пространстве
struct point2D
{
   double x, y;

   point2D(double _x = 0, double _y = 0) :
      x(_x), y(_y)
   {
   };
};

// Описание подобласти (параметры уравнения)
struct area
{
   double x_start, x_end;
   double y_start, y_end;

   double lambda, sigma, hi;
};

// Опорный элемент (прямоугольник) в 2D
struct omega
{
   point2D start_point; // (xk,  yk)
   point2D end_point;   // (xk1, yk1)
};

#endif