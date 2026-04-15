// ================  TESTS.CPP ================
#include "tests.hpp"

// Индекс активного набора тест-функций (0..3)
uint8_t ifu = 3;

/*
  Гармоническая задача (вар. 5): 2D, декартовы координаты, билинейные базис. функции.
  Уравнение:
      -div(λ grad u_s) - ω²χ u_s - ωσ u_c = f_s
      -div(λ grad u_c) - ω²χ u_c + ωσ u_s = f_c

  Тест 0: u_s = x+y,     u_c = x-y       (линейные, Δu = 0)
  Тест 1: u_s = x²+y²,   u_c = x²-y²    (квадратичные)
  Тест 2: u_s = x³+y³,   u_c = x³-y³    (кубические)
  Тест 3: u_s = sin(x+y),u_c = exp(x+y) (нелинейные)
*/
tests::tests()
{
   funcs.resize(4);

   // ---- Тест 0: полиномы 1-й степени (Δu = 0) ----
   funcs[0].us = [](double x, double y) { return x + y; };
   funcs[0].uc = [](double x, double y) { return x - y; };
   funcs[0].fs = [this](double x, double y, double lambda, double sigma, double hi, double w)
      {
         return -w * sigma * funcs[0].uc(x, y)
            - w * w * hi * funcs[0].us(x, y);
      };
   funcs[0].fc = [this](double x, double y, double lambda, double sigma, double hi, double w)
      {
         return w * sigma * funcs[0].us(x, y)
            - w * w * hi * funcs[0].uc(x, y);
      };

   // ---- Тест 1: полиномы 2-й степени ----
   // Δ(x²+y²) = 4,   Δ(x²-y²) = 0
   funcs[1].us = [](double x, double y) { return x * x + y * y; };
   funcs[1].uc = [](double x, double y) { return x * x - y * y; };
   funcs[1].fs = [this](double x, double y, double lambda, double sigma, double hi, double w)
      {
         return -4.0 * lambda
            - w * sigma * funcs[1].uc(x, y)
            - w * w * hi * funcs[1].us(x, y);
      };
   funcs[1].fc = [this](double x, double y, double lambda, double sigma, double hi, double w)
      {
         return  w * sigma * funcs[1].us(x, y)
            - w * w * hi * funcs[1].uc(x, y);
      };

   // ---- Тест 2: полиномы 3-й степени ----
   // Δ(x³+y³) = 6x+6y,  Δ(x³-y³) = 6x-6y
   funcs[2].us = [](double x, double y) { return x * x * x + y * y * y; };
   funcs[2].uc = [](double x, double y) { return x * x * x - y * y * y; };
   funcs[2].fs = [this](double x, double y, double lambda, double sigma, double hi, double w)
      {
         return -lambda * 6.0 * (x + y)
            - w * sigma * funcs[2].uc(x, y)
            - w * w * hi * funcs[2].us(x, y);
      };
   funcs[2].fc = [this](double x, double y, double lambda, double sigma, double hi, double w)
      {
         return -lambda * 6.0 * (x - y)
            + w * sigma * funcs[2].us(x, y)
            - w * w * hi * funcs[2].uc(x, y);
      };

   // ---- Тест 3: нелинейные функции ----
   // Δ(sin(x+y)) = -2 sin(x+y),  Δ(exp(x+y)) = 2 exp(x+y)
   funcs[3].us = [](double x, double y) { return sin(x + y); };
   funcs[3].uc = [](double x, double y) { return exp(x + y); };
   funcs[3].fs = [this](double x, double y, double lambda, double sigma, double hi, double w)
      {
         return 2.0 * lambda * sin(x + y)
            - w * sigma * funcs[3].uc(x, y)
            - w * w * hi * funcs[3].us(x, y);
      };
   funcs[3].fc = [this](double x, double y, double lambda, double sigma, double hi, double w)
      {
         return -2.0 * lambda * exp(x + y)
            + w * sigma * funcs[3].us(x, y)
            - w * w * hi * funcs[3].uc(x, y);
      };
}

// ----------------------------------------------------------------
void tests::calc_exact(space_grid &grid)
{
   exact.resize(2 * grid.get_nodes_count());

   for (uint32_t i = 0; i < grid.get_nodes_count(); i++)
   {
      point2D pt = grid.get_point(i);
      exact[2 * i] = funcs[ifu].us(pt.x, pt.y);
      exact[2 * i + 1] = funcs[ifu].uc(pt.x, pt.y);
   }
}

// ----------------------------------------------------------------
// Заголовок таблицы (сравниваем LU и LOS_LU)
static void write_header(std::ofstream &out, const std::string &param)
{
   out << std::left
      << std::setw(12) << param
      << std::setw(12) << "|LU_time"
      << std::setw(12) << "|LOS_time"
      << std::setw(14) << "|LOS_iters"
      << std::setw(18) << "|LOS_residual"
      << std::setw(14) << "|LU_error"
      << std::setw(14) << "|LOS_error"
      << "\n";
   out << std::string(96, '-') << "\n";
}

// ----------------------------------------------------------------
void tests::omega_tests()
{
   space_grid_generator sgg;
   space_grid *sg = nullptr;
   sgg.build_mesh(sg, funcs[ifu].us, funcs[ifu].uc);

   calc_exact(*sg);
   double exact_norm = norm(exact);

   double omegas[] = { 1e-4, 1e-3, 1e-2, 1e-1, 1, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9 };

   std::ofstream out(directory + "tests\\omega_tests.txt");
   write_header(out, "omega");

   std::cout << "=== omega tests ===\n";

   for (const auto &w : omegas)
   {
      mfe solver(*sg);
      solver.set_w(w);
      solver.assembly_global_matrix_and_vector(funcs[ifu].fs, funcs[ifu].fc);

      result r1 = solver.solve(method::LU);
      double err_lu = norm(solver.q - exact) / exact_norm;

      result r2 = solver.solve(method::LOS_LU);
      double err_los = norm(solver.q - exact) / exact_norm;

      out << std::left
         << std::setw(12) << w
         << "|" << std::setw(11) << r1.time.count()
         << "|" << std::setw(11) << r2.time.count()
         << "|" << std::setw(13) << r2.iters
         << "|" << std::setw(17) << r2.residual
         << "|" << std::setw(13) << err_lu
         << "|" << std::setw(13) << err_los
         << "\n";

      std::cout << "w=" << w << "  LU_err=" << err_lu << "  LOS_err=" << err_los << "\n";
   }
   std::cout << "\n";

   delete sg;
}

// ----------------------------------------------------------------
void tests::lambda_tests()
{
   space_grid_generator sgg;
   space_grid *sg = nullptr;
   sgg.build_mesh(sg, funcs[ifu].us, funcs[ifu].uc);

   calc_exact(*sg);
   double exact_norm = norm(exact);

   double lambdas[] = { 1e2, 1e3, 1e4, 1e5, 8e5 };

   std::ofstream out(directory + "tests\\lambda_tests.txt");
   write_header(out, "lambda");

   std::cout << "=== lambda tests ===\n";

   for (const auto &lam : lambdas)
   {
      sg->set_lambda(lam);

      mfe solver(*sg);
      solver.set_w(10);
      solver.assembly_global_matrix_and_vector(funcs[ifu].fs, funcs[ifu].fc);

      result r1 = solver.solve(method::LU);
      double err_lu = norm(solver.q - exact) / exact_norm;

      result r2 = solver.solve(method::LOS_LU);
      double err_los = norm(solver.q - exact) / exact_norm;

      out << std::left
         << std::setw(12) << lam
         << "|" << std::setw(11) << r1.time.count()
         << "|" << std::setw(11) << r2.time.count()
         << "|" << std::setw(13) << r2.iters
         << "|" << std::setw(17) << r2.residual
         << "|" << std::setw(13) << err_lu
         << "|" << std::setw(13) << err_los
         << "\n";

      std::cout << "lambda=" << lam << "  LU_err=" << err_lu << "  LOS_err=" << err_los << "\n";
   }
   std::cout << "\n";

   delete sg;
}

// ----------------------------------------------------------------
void tests::sigma_tests()
{
   space_grid_generator sgg;
   space_grid *sg = nullptr;
   sgg.build_mesh(sg, funcs[ifu].us, funcs[ifu].uc);

   calc_exact(*sg);
   double exact_norm = norm(exact);

   double sigmas[] = { 0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8 };

   std::ofstream out(directory + "tests\\sigma_tests.txt");
   write_header(out, "sigma");

   std::cout << "=== sigma tests ===\n";

   for (const auto &sig : sigmas)
   {
      sg->set_sigma(sig);

      mfe solver(*sg);
      solver.set_w(10);
      solver.assembly_global_matrix_and_vector(funcs[ifu].fs, funcs[ifu].fc);

      result r1 = solver.solve(method::LU);
      double err_lu = norm(solver.q - exact) / exact_norm;

      result r2 = solver.solve(method::LOS_LU);
      double err_los = norm(solver.q - exact) / exact_norm;

      out << std::left
         << std::setw(12) << sig
         << "|" << std::setw(11) << r1.time.count()
         << "|" << std::setw(11) << r2.time.count()
         << "|" << std::setw(13) << r2.iters
         << "|" << std::setw(17) << r2.residual
         << "|" << std::setw(13) << err_lu
         << "|" << std::setw(13) << err_los
         << "\n";

      std::cout << "sigma=" << sig << "  LU_err=" << err_lu << "  LOS_err=" << err_los << "\n";
   }
   std::cout << "\n";

   delete sg;
}

// ----------------------------------------------------------------
void tests::hi_tests()
{
   space_grid_generator sgg;
   space_grid *sg = nullptr;
   sgg.build_mesh(sg, funcs[ifu].us, funcs[ifu].uc);

   calc_exact(*sg);
   double exact_norm = norm(exact);

   double his[] = { 8.81e-12, 1e-12, 1e-11, 1e-10 };

   std::ofstream out(directory + "tests\\hi_tests.txt");
   write_header(out, "hi");

   std::cout << "=== hi tests ===\n";

   for (const auto &hi : his)
   {
      sg->set_hi(hi); // исправлен баг: было set_sigma(hi)

      mfe solver(*sg);
      solver.set_w(10);
      solver.assembly_global_matrix_and_vector(funcs[ifu].fs, funcs[ifu].fc);

      result r1 = solver.solve(method::LU);
      double err_lu = norm(solver.q - exact) / exact_norm;

      result r2 = solver.solve(method::LOS_LU);
      double err_los = norm(solver.q - exact) / exact_norm;

      out << std::left
         << std::setw(12) << hi
         << "|" << std::setw(11) << r1.time.count()
         << "|" << std::setw(11) << r2.time.count()
         << "|" << std::setw(13) << r2.iters
         << "|" << std::setw(17) << r2.residual
         << "|" << std::setw(13) << err_lu
         << "|" << std::setw(13) << err_los
         << "\n";

      std::cout << "hi=" << hi << "  LU_err=" << err_lu << "  LOS_err=" << err_los << "\n";
   }
   std::cout << "\n";

   delete sg;
}

// ----------------------------------------------------------------
void tests::run()
{
   omega_tests();
   lambda_tests();
   sigma_tests();
   hi_tests();
}