#include "testing_module.h"

void testing_module::set_functions()
{
#pragma region Выделяем память
   u.resize(6);
   lambda.resize(5);
   f.resize(6);

   for (uint32_t i = 0; i < f.size(); i++)
      f[i].resize(5);

#pragma endregion

#pragma region Неизвестная функция

   u[0] = [](double x, double t) { return x + t; };
   u[1] = [](double x, double t) { return x * x + t; };
   u[2] = [](double x, double t) { return t * t; };
   u[3] = [](double x, double t) { return x * x * x + t; };
   u[4] = [](double x, double t) { return x * x * x * x + t; };
   u[5] = [](double x, double t) { return exp(x) + sin(t); };

   u_names = {
       "x+t",
       "x^2+t",
       "x+t^2",
       "x^3+t",
       "x^4+t",
       "e^x+sin(t)"
   };
#pragma endregion

#pragma region Функция лямбда

   lambda[0] = [](double q2, double q1, double x) { return 1; };
   lambda[1] = [](double q2, double q1, double x) { return x + 1; };
   lambda[2] = [](double q2, double q1, double x) { return q2 - q1 + 1; };
   lambda[3] = [](double q2, double q1, double x) { return (q2 - q1) * (q2 - q1) + 4; };
   lambda[4] = [](double q2, double q1, double x) { return exp(q2 - q1); };

   lambda_names = {
       "1",
       "x+1",
       "du/dx+1",
       "(du/dx)^2+4",
       "e^(du/dx)"
   };

#pragma endregion

#pragma region Функция правой части

   f[0][0] = [](double x, double t, double sigma) { return sigma; };
   f[0][1] = [](double x, double t, double sigma) { return sigma - 1; };
   f[0][2] = [](double x, double t, double sigma) { return sigma; };
   f[0][3] = [](double x, double t, double sigma) { return sigma; };
   f[0][4] = [](double x, double t, double sigma) { return sigma; };

   f[1][0] = [](double x, double t, double sigma) { return sigma - 2; };
   f[1][1] = [](double x, double t, double sigma) { return sigma - 2 - 4 * x; };
   f[1][2] = [](double x, double t, double sigma) { return sigma - 2 - 8 * x; };
   f[1][3] = [](double x, double t, double sigma) { return sigma - 8 - 24 * x * x; };
   f[1][4] = [](double x, double t, double sigma) { return sigma - 2 * exp(2 * x) * (2 * x + 1); };

   f[2][0] = [](double x, double t, double sigma) { return 2 * t * sigma; };
   f[2][1] = [](double x, double t, double sigma) { return 2 * t * sigma - 1; };
   f[2][2] = [](double x, double t, double sigma) { return 2 * t * sigma; };
   f[2][3] = [](double x, double t, double sigma) { return 2 * t * sigma; };
   f[2][4] = [](double x, double t, double sigma) { return 2 * t * sigma; };

   f[3][0] = [](double x, double t, double sigma) { return sigma - 6 * x; };
   f[3][1] = [](double x, double t, double sigma) { return -9 * x * x - 6 * x + sigma; };
   f[3][2] = [](double x, double t, double sigma) { return sigma - 36 * x * x * x - 6 * x; };
   f[3][3] = [](double x, double t, double sigma) { return sigma - 162 * x * x * x * x * x - 24 * x; };
   f[3][4] = [](double x, double t, double sigma) { return sigma - exp(3 * x * x) * 6 * x * (3 * x * x + 1); };

   f[4][0] = [](double x, double t, double sigma) { return sigma - 12 * x * x; };
   f[4][1] = [](double x, double t, double sigma) { return sigma - 16 * x * x * x - 12 * x * x; };
   f[4][2] = [](double x, double t, double sigma) { return sigma - 96 * pow(x, 5) - 12 * x * x; };
   f[4][3] = [](double x, double t, double sigma) { return sigma - 576 * pow(x, 8) - 48 * x * x; };
   f[4][4] = [](double x, double t, double sigma) { return sigma - exp(4 * x * x * x) * 12 * x * x * (4 * x * x * x + 1); };

   f[5][0] = [](double x, double t, double sigma) { return -exp(x) + sigma * cos(t); };
   f[5][1] = [](double x, double t, double sigma) { return -exp(x) * (2 + x) + cos(t) * sigma; };
   f[5][2] = [](double x, double t, double sigma) { return -exp(x) * (2 * exp(x) + 1) + cos(t) * sigma; };
   f[5][3] = [](double x, double t, double sigma) { return -exp(x) * (3 * exp(2 * x) + 4) + cos(t) * sigma; };
   f[5][4] = [](double x, double t, double sigma) { return -exp(exp(x) + x) * (exp(x) + 1) + cos(t) * sigma; };

#pragma endregion
}

void testing_module::save(std::ofstream &file, std::pair<uint32_t, double> &res)
{
   if (file.is_open())
   {
      file << std::left << std::scientific
         << std::setw(15) << "Iterations: "
         << res.first
         << std::endl
         << std::setw(15) << "Error: "
         << res.second
         << std::endl;
      file << "----------------------" << std::endl;
   }
   else
      std::cerr << "File is not open" << std::endl;
}

void testing_module::run()
{
   mesh mesh;
   mesh_generator mg;

   mg.build_mesh(mesh);

   mfe mfe(mesh);

   //------------------------------------------------------------------------------
   std::ofstream table("table.txt", std::ios::app);
   table.precision(4);

   table << std::left << std::setw(15) << "u/lambda" << "|";
   for (uint32_t i = 0; i < lambda_names.size(); i++)
      table << std::setw(20) << lambda_names[i];
   table << std::endl;
   table << "--------------------------------------------------------";
   table << "--------------------------------------------------------" << std::endl;
   //------------------------------------------------------------------------------

   for (uint32_t i = 0; i < u.size(); i++)
   {
      //------------------------------------------------------------------------------
      table << std::left << std::setw(15) << u_names[i] << "|";
      //------------------------------------------------------------------------------

      std::ofstream tres(tests_directory + "test_" + std::to_string(i) + ".txt", std::ios::out | std::ios::trunc);

      for (uint32_t j = 0; j < lambda.size(); j++)
      {
         mfe.set_functions(u[i], f[i][j], lambda[j]);

         auto res = mfe.solve(mesh, 5000, 1e-13);

         std::stringstream ss;
         std::string str;

         ss << std::scientific << res.second;
         ss >> str;
         str += "(" + std::to_string(res.first) + ")";

         //-----------------------------------------------------------------------------
         table << std::scientific << std::left << std::setw(20) << str;
         //-----------------------------------------------------------------------------

         save(tres, res);
      }
      table << std::endl;

      tres.close();
   }

   table << std::endl << std::endl;
   table.close();

   u.clear();
   lambda.clear();
   f.clear();
}