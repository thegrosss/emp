// ================  TESTS.CPP ================
#include "tests.hpp"

uint8_t ifu = 3;

tests::tests()
{
   funcs.resize(4);

   funcs[0].us = [](double x, double y, double z) { return x + y + z; };
   funcs[0].uc = [](double x, double y, double z) { return x - y - z; };
   funcs[0].fs = [this](double x, double y, double z, double lambda, double sigma, double hi, double w)
      {
         return -w * sigma * funcs[0].uc(x, y, z) - w * w * hi * funcs[0].us(x, y, z);
      };
   funcs[0].fc = [this](double x, double y, double z, double lambda, double sigma, double hi, double w)
      {
         return w * sigma * funcs[0].us(x, y, z) - w * w * hi * funcs[0].uc(x, y, z);
      };

   funcs[1].us = [](double x, double y, double z) { return x * x + y * y + z * z; };
   funcs[1].uc = [](double x, double y, double z) { return 3 * x * x - 2 * y * y + z * z; };
   funcs[1].fs = [this](double x, double y, double z, double lambda, double sigma, double hi, double w)
      {
         return -6.0 * lambda - w * sigma * funcs[1].uc(x, y, z) - w * w * hi * funcs[1].us(x, y, z);
      };
   funcs[1].fc = [this](double x, double y, double z, double lambda, double sigma, double hi, double w)
      {
         return -4.0 * lambda + w * sigma * funcs[1].us(x, y, z) - w * w * hi * funcs[1].uc(x, y, z);
      };

   funcs[2].us = [](double x, double y, double z) { return x * x * x + y * y * y + z * z * z; };
   funcs[2].uc = [](double x, double y, double z) { return 2 * x * x * x - y * y * y + 3 * z * z * z; };
   funcs[2].fs = [this](double x, double y, double z, double lambda, double sigma, double hi, double w)
      {
         return -lambda * 6 * (x + y + z) - w * sigma * funcs[2].uc(x, y, z) - w * w * hi * funcs[2].us(x, y, z);
      };
   funcs[2].fc = [this](double x, double y, double z, double lambda, double sigma, double hi, double w)
      {
         return -lambda * 6 * (2 * x - y + 3 * z) + w * sigma * funcs[2].us(x, y, z) - w * w * hi * funcs[2].uc(x, y, z);
      };

   funcs[3].us = [](double x, double y, double z) { return sin(x + y + z); };
   funcs[3].uc = [](double x, double y, double z) { return exp(x + y) - z; };
   funcs[3].fs = [this](double x, double y, double z, double lambda, double sigma, double hi, double w)
      {
         return lambda * sin(x + y + z) - w * sigma * funcs[3].uc(x, y, z) - w * w * hi * funcs[3].us(x, y, z);
      };
   funcs[3].fc = [this](double x, double y, double z, double lambda, double sigma, double hi, double w)
      {
         return -2.0 * lambda * exp(x + y) + w * sigma * funcs[3].us(x, y, z) - w * w * hi * funcs[3].uc(x, y, z);
      };
}


void tests::calc_exact(space_grid &grid)
{
   exact.resize(2 * grid.get_nodes_count());

   for (uint32_t i = 0; i < grid.get_nodes_count(); i++)
   {
      point3D point = grid.get_point(i);

      exact[2 * i] = funcs[ifu].us(point.x, point.y, point.z);
      exact[2 * i + 1] = funcs[ifu].uc(point.x, point.y, point.z);
   }
}


void tests::omega_tests()
{
   space_grid_generator sgg;
   space_grid *sg = nullptr;

   sgg.build_mesh(sg, funcs[ifu].us, funcs[ifu].uc);

   calc_exact(*sg);
   double exact_norm = norm(exact);

   double omegas[] = { 1e-4, 1e-3, 1e-2, 1e-1, 1, 1e1, 1e2, 1e-3, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9 };

   result res1, res2, res3;

   std::ofstream wres(directory + "tests\\" + "omega_tests.txt");

   wres << std::left << std::setw(10) << "omega"
      << std::setw(10) << "|LU_time" << std::setw(10) << "|LOS_time" << std::setw(14) << "|BCGSTAB_time"
      << std::setw(14) << "|LOS_iters" << std::setw(16) << "|BCGSTAB_iters"
      << std::setw(16) << "|LOS_residual" << std::setw(16) << "|BSG_residual"
      << std::setw(16) << "|LOS_error" << std::setw(16) << "|BCG_error" << std::endl;
   wres << "--------------------------------------------------------------------";
   wres << "------------------------------------------------------------------" << std::endl;

   std::cout << "omega" << std::endl;

   for (const auto& w : omegas)
   {
      //double w = 1;

      mfe mfe(*sg);
      mfe.set_w(w);
      mfe.assembly_global_matrix_and_vector(funcs[ifu].fs, funcs[ifu].fc);

      res1 = mfe.solve(method::LU);
      auto dif = mfe.q - exact;
      double err0 = norm(dif) / exact_norm;
      std::cout << err0 << std::endl;

      res2 = mfe.solve(method::LOS_LU);
      dif = mfe.q - exact;
      double err1 = norm(dif) / exact_norm;

      res3 = mfe.solve(method::BCGSTAB_LU);
      dif = mfe.q - exact;
      double err2 = norm(dif) / exact_norm;

      wres << std::left << std::setw(10) << w
         << "|" << std::setw(9) << res1.time.count() << "|" << std::setw(9) << res2.time.count() << "|" << std::setw(13) << res3.time.count()
         << "|" << std::setw(13) << res2.iters << "|" << std::setw(15) << res3.iters
         << "|" << std::setw(15) << res2.residual << "|" << std::setw(15) << res3.residual
         << "|" << std::setw(15) << err1 << "|" << std::setw(15) << err2 << std::endl;

      mfe.~mfe();
   }
   std::cout << std::endl;
}

void tests::lambda_tests()
{
   space_grid_generator sgg;
   space_grid *sg = nullptr;

   sgg.build_mesh(sg, funcs[ifu].us, funcs[ifu].uc);

   calc_exact(*sg);
   double exact_norm = norm(exact);

   double lambdas[] = { 1e2, 1e3, 1e4, 1e5, 8e5 };

   result res1, res2, res3;

   std::ofstream wres(directory + "tests\\" + "lambda_tests.txt");
   std::ofstream wresdif(directory + "tests\\" + "lambda_tests_diffs.txt");

   wres << std::left << std::setw(10) << "lambda"
      << std::setw(10) << "|LU_time" << std::setw(10) << "|LOS_time" << std::setw(14) << "|BCGSTAB_time"
      << std::setw(14) << "|LOS_iters" << std::setw(16) << "|BCGSTAB_iters"
      << std::setw(16) << "|LOS_residual" << std::setw(16) << "|BSG_residual"
      << std::setw(16) << "|LOS_error" << std::setw(16) << "|BCG_error" << std::endl;
   wres << "--------------------------------------------------------------------";
   wres << "------------------------------------------------------------------" << std::endl;

   std::cout << "lambda" << std::endl;

   for (const auto &lambda : lambdas)
   {
      sg->set_lambda(lambda);

      mfe mfe(*sg);
      mfe.set_w(10);
      mfe.assembly_global_matrix_and_vector(funcs[ifu].fs, funcs[ifu].fc);

      res1 = mfe.solve(method::LU);
      auto dif = mfe.q - exact;
      double err0 = norm(dif) / exact_norm;
      std::cout << err0 << std::endl;

      res2 = mfe.solve(method::LOS_LU);
      dif = mfe.q - exact;
      double err1 = norm(dif) / exact_norm;

      res3 = mfe.solve(method::BCGSTAB_LU);
      dif = mfe.q - exact;
      double err2 = norm(dif) / exact_norm;


      wres << std::left << std::setw(10) << lambda
         << "|" << std::setw(9) << res1.time.count() << "|" << std::setw(9) << res2.time.count() << "|" << std::setw(13) << res3.time.count()
         << "|" << std::setw(13) << res2.iters << "|" << std::setw(15) << res3.iters <<
         "|" << std::setw(15) << res2.residual << "|" << std::setw(15) << res3.residual
         << "|" << std::setw(15) << err1 << "|" << std::setw(15) << err2 << std::endl;

      mfe.~mfe();
   }
   std::cout << std::endl;
}

void tests::sigma_tests()
{
   space_grid_generator sgg;
   space_grid *sg = nullptr;

   sgg.build_mesh(sg, funcs[ifu].us, funcs[ifu].uc);

   calc_exact(*sg);
   double exact_norm = norm(exact);

   double sigmas[] = { 0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8 };

   result res1, res2, res3;

   std::ofstream wres(directory + "tests\\" + "sigma_tests.txt");
   std::ofstream wresdif(directory + "tests\\" + "sigma_tests_diffs.txt");

   wres << std::left << std::setw(10) << "sigma"
      << std::setw(10) << "|LU_time" << std::setw(10) << "|LOS_time" << std::setw(14) << "|BCGSTAB_time"
      << std::setw(14) << "|LOS_iters" << std::setw(16) << "|BCGSTAB_iters"
      << std::setw(16) << "|LOS_residual" << std::setw(16) << "|BSG_residual"
      << std::setw(16) << "|LOS_error" << std::setw(16) << "|BCG_error" << std::endl;
   wres << "--------------------------------------------------------------------";
   wres << "------------------------------------------------------------------" << std::endl;

   std::cout << "sigma" << std::endl;

   for (const auto &sigma : sigmas)
   {
      sg->set_sigma(sigma);

      mfe mfe(*sg);
      mfe.set_w(10);
      mfe.assembly_global_matrix_and_vector(funcs[ifu].fs, funcs[ifu].fc);

      res1 = mfe.solve(method::LU);
      auto dif = mfe.q - exact;
      double err0 = norm(dif) / exact_norm;
      std::cout << err0 << std::endl;

      res2 = mfe.solve(method::LOS_LU);
      dif = mfe.q - exact;
      double err1 = norm(dif) / exact_norm;

      res3 = mfe.solve(method::BCGSTAB_LU);
      dif = mfe.q - exact;
      double err2 = norm(dif) / exact_norm;

      wres << std::left << std::setw(10) << sigma
         << "|" << std::setw(9) << res1.time.count() << "|" << std::setw(9) << res2.time.count() << "|" << std::setw(13) << res3.time.count()
         << "|" << std::setw(13) << res2.iters << "|" << std::setw(15) << res3.iters <<
         "|" << std::setw(15) << res2.residual << "|" << std::setw(15) << res3.residual
         << "|" << std::setw(15) << err1 << "|" << std::setw(15) << err2 << std::endl;

      mfe.~mfe();
   }
   std::cout << std::endl;
}

void tests::hi_tests()
{
   space_grid_generator sgg;
   space_grid *sg = nullptr;

   sgg.build_mesh(sg, funcs[ifu].us, funcs[ifu].uc);

   calc_exact(*sg);
   double exact_norm = norm(exact);

   double his[] = { 8.81e-12, 1e-12, 1e-11, 1e-10 };

   result res1, res2, res3;

   std::ofstream wres(directory + "tests\\" + "hi_tests.txt");
   std::ofstream wresdif(directory + "tests\\" + "hi_tests_diffs.txt");

   wres << std::left << std::setw(10) << "hi"
      << std::setw(10) << "|LU_time" << std::setw(10) << "|LOS_time" << std::setw(14) << "|BCGSTAB_time"
      << std::setw(14) << "|LOS_iters" << std::setw(16) << "|BCGSTAB_iters"
      << std::setw(16) << "|LOS_residual" << std::setw(16) << "|BSG_residual"
      << std::setw(16) << "|LOS_error" << std::setw(16) << "|BCG_error" << std::endl;
   wres << "--------------------------------------------------------------------";
   wres << "------------------------------------------------------------------" << std::endl;

   std::cout << "hi" << std::endl;

   for (const auto &hi : his)
   {
      sg->set_hi(hi);

      mfe mfe(*sg);
      mfe.set_w(10);
      mfe.assembly_global_matrix_and_vector(funcs[ifu].fs, funcs[ifu].fc);

      res1 = mfe.solve(method::LU);
      auto dif = mfe.q - exact;
      double err0 = norm(dif) / exact_norm;
      std::cout << err0 << std::endl;

      res2 = mfe.solve(method::LOS_LU);
      dif = mfe.q - exact;
      double err1 = norm(dif) / exact_norm;

      res3 = mfe.solve(method::BCGSTAB_LU);
      dif = mfe.q - exact;
      double err2 = norm(dif) / exact_norm;

      wres << std::left << std::setw(10) << hi
         << "|" << std::setw(9) << res1.time.count() << "|" << std::setw(9) << res2.time.count() << "|" << std::setw(13) << res3.time.count()
         << "|" << std::setw(13) << res2.iters << "|" << std::setw(15) << res3.iters <<
         "|" << std::setw(15) << res2.residual << "|" << std::setw(15) << res3.residual
         << "|" << std::setw(15) << err1 << "|" << std::setw(15) << err2 << std::endl;

      mfe.~mfe();
   }

   std::cout << std::endl;
}


void tests::run()
{
   omega_tests();

   lambda_tests();

   sigma_tests();

   hi_tests();
}
