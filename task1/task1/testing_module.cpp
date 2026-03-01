#include "testing_module.h"

// Задаем тестовые функции
void testing_module::set_functions()
{
   functions.resize(5);

   // 1-ый элемент пары - неизввестная функция u
   // 2-ой функция правой части f
   functions[0] = std::make_pair
   (
      [](double x, double y) { return x + 3 * y; },
      [](double x, double y, double c1, double c2) { return c2 * (x + 3 * y); }
   );

   functions[1] = std::make_pair
   (
      [](double x, double y) { return 2 * x * x - 3 * y * y; },
      [](double x, double y, double c1, double c2) { return 2 * c1 + c2 * (2 * x * x - 3 * y * y); }
   );

   functions[2] = std::make_pair
   (
      [](double x, double y) { return x * x * x + y * y * y; },
      [](double x, double y, double c1, double c2) { return  c1 * (-6 * x - 6 * y) + c2 * (x * x * x + y * y * y); }
   );

   functions[3] = std::make_pair
   (
      [](double x, double y) { return x * x * x * x - 10 * y * y * y * y; },
      [](double x, double y, double c1, double c2) { return c1 * (-12 * x * x + 120 * y * y) + c2 * (x * x * x * x - 10 * y * y * y * y); }
   );

   functions[4] = std::make_pair
   (
      [](double x, double y) { return exp(x) + y; },
      [](double x, double y, double c1, double c2) { return -exp(x) * c1 + c2 * (exp(x) + y); }
   );
}

// Запускаем все тесты
void testing_module::run_tests()
{
   mesh_generator mg;
   mesh mesh;
   fdm fdm;
   std::vector<double> u;

   mg.build_mesh(mesh);

   for (uint32_t test = 0; test < functions.size(); test++)
   {
      calc_exact(functions[test].first, mesh);
      fdm.exact = exact;

      fdm.mesh_to_slae(mesh, functions[test].first, functions[test].second);

      auto result = fdm.calculate(u);

      save(test + 1, result, u);
   }

   mesh.clear();
}

// Сохраняем результат в файл
void testing_module::save(uint32_t test_num,
   std::pair<uint32_t, double> &result,
   const std::vector<double> &u)
{
   std::ofstream out("result/test" + std::to_string(test_num) + ".txt");

   out << std::left
      << std::setw(16) << "iterations: " << result.first << std::endl
      << std::setw(16) << "residual: " << result.second << std::endl
      << std::setw(16) << "relative error: " << error(u) << std::endl << std::endl;

   out.precision(13);
   out << std::left << std::setw(20) << "u"
      << std::setw(20) << "u*"
      << std::setw(20) << "|u-u*|" << std::endl;

   for (uint32_t i = 0; i < u.size(); i++)
   {
      out << std::left << std::setw(20) << u[i]
         << std::setw(20) << exact[i]
         << std::setw(20) << abs(exact[i] - u[i]) << std::endl;
   }
   out.close();
}

// Считаем относительную погрешность
double testing_module::error(const std::vector<double> &q)
{
   double sum_sqr_diff = 0.0;
   double sum_sqr = 0.0;
   for (uint32_t i = 0; i < q.size(); i++)
   {
      sum_sqr_diff += (exact[i] - q[i]) * (exact[i] - q[i]);
      sum_sqr += exact[i] * exact[i];
   }
   return sqrt(sum_sqr_diff) / sqrt(sum_sqr);
}

// Считаем точное значение функции
void testing_module::calc_exact(const func2D_u &u, mesh &mesh)
{
   exact.resize(mesh.size());

   for (uint32_t i = 0; i < exact.size(); i++)
   {
      if (mesh[i].type != -1)
         exact[i] = u(mesh[i].p.x, mesh[i].p.y);
      else
         exact[i] = 0.0;
   }
}