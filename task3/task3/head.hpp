// ================ HEAD.HPP ================
#ifndef HEAD_HPP
#define HEAD_HPP

#include <iostream>
#include <vector>
#include <array>
#include <set>
#include <fstream>
#include <iomanip>
#include <functional>
#include <string>
#include <chrono>
#include "json.hpp"

const std::string directory = "C:\\Users\\kiril\\Desktop\\НГТУ\\3 курс\\6 семестр\\Уравнения математической физики\\task3\\task3\\";

typedef std::vector<double> dvector;

// Функция двух пространственных переменных
typedef std::function<double(double, double)> function2D;

// Правая часть / граничное условие: (x, y, lambda, sigma, hi, w)
typedef std::function<double(double, double, double, double, double, double)> r_function2D;

#endif