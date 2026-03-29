#pragma once
#include <iostream>
#include <vector>
#include <functional>
#include <fstream>
#include <iomanip>
#include <string>
#include <sstream>

typedef std::function<double(double, double)> function;
typedef std::function<double(double, double, double)> function_f;
typedef std::function<double(double, double, double)> function_lambda;
typedef std::function<double(double)> min_func;

const std::string directory = "C:\\Users\\kiril\\Desktop\\НГТУ\\3 курс\\6 семестр\\Уравнения математической физики\\Лаба 2\\task2\\task2\\";
const std::string tests_directory = directory + "tests\\";

// Посчитать норму вектора
inline double norm(const std::vector<double> &v)
{
   double scalar = 0.0;

   for (const auto &it : v)
      scalar += it * it;

   return sqrt(scalar);
}