#include "mfe.h"
#include "json_reader.h"

mfe::mfe(mesh &mesh)
{
	JsonConfig config = JsonConfig::from_file(directory + "data.json");

   global_A = new matrix(mesh.get_fe_count() + 1, 3);
   A = nullptr;
   global_f.resize(mesh.get_fe_count() + 1);
   q.resize(mesh.get_fe_count() + 1);
   q_prev.resize(mesh.get_fe_count() + 1);

   if (config.mfe_method == "SIMPLE_ITERATION")
      method_to_solve = method::SIMPLE_ITERATION;
   else if (config.mfe_method == "NEWTON")
      method_to_solve = method::NEWTON;
   else
		throw "Unknown method to solve";

   use_relax = config.use_relax;

   local_f.resize(2);
   local_A.resize(2);
   for (auto &row : local_A)
      row.resize(2);
}

void mfe::set_functions(const function &_u, const function_f &_f, const function_lambda &_lambda)
{
   u = _u;
   f = _f;
   lambda = _lambda;
}

// Сборка локальной матрицы
void mfe::build_local_matrix(const finite_elem &elem, const double delta_t)
{
   double h = elem.x_next - elem.x;

   double lambda1 = lambda(q[elem.number + 1] / h, q[elem.number] / h, elem.x);
   double lambda2 = lambda(q[elem.number + 1] / h, q[elem.number] / h, elem.x_next);

   // Матрица жесткости
   double coef = (lambda1 + lambda2) / (2.0 * h);
   local_A[0][0] = local_A[1][1] = coef;
   local_A[0][1] = local_A[1][0] = -coef;

   // Матрица жесткости + матрица масс
   coef = elem.sigma * h / (6.0 * delta_t);
   local_A[0][0] += 2 * coef;
   local_A[1][1] += 2 * coef;
   local_A[0][1] += coef;
   local_A[1][0] += coef;
}

// Сборка локального вектора
void mfe::build_local_vector(const finite_elem &elem, const double t, const double delta_t)
{
   double h = elem.x_next - elem.x;

   double f1 = f(elem.x, t, elem.sigma);
   double f2 = f(elem.x_next, t, elem.sigma);

   local_f[0] = h * (2 * f1 + f2) / 6.0 + elem.sigma * h * (2 * q_prev[elem.number] + q_prev[elem.number + 1]) / (6.0 * delta_t);
   local_f[1] = h * (f1 + 2 * f2) / 6.0 + elem.sigma * h * (q_prev[elem.number] + 2 * q_prev[elem.number + 1]) / (6.0 * delta_t);
}

// Сборка глоабальной СЛАУ (глобального вектора и матрицы)
void mfe::assembly_global_slae(mesh &mesh, const double t, const double delta_t)
{
   global_A->set_elem_to_null();
   global_f.clear();
   global_f.resize(q.size());

   if (method_to_solve == method::NEWTON)
   {
      A->set_elem_to_null();
      nl_f.clear();
      nl_f.resize(q.size());
   }

   for (uint32_t ielem = 0; ielem < mesh.get_fe_count(); ielem++)
   {
      // Собираем локальную матрицу
      build_local_matrix(mesh[ielem], delta_t);

      // Собираем локальынй вектор
      build_local_vector(mesh[ielem], t, delta_t);

      if (method_to_solve == method::NEWTON)
      {
         A->add_to_elem(0, ielem, local_A[0][0]);
         A->add_to_elem(0, ielem + 1, local_A[1][1]);
         A->add_to_elem(1, ielem, local_A[1][0]);
         A->add_to_elem(2, ielem, local_A[0][1]);

         nl_f[ielem] += local_f[0];
         nl_f[ielem + 1] += local_f[1];

         linearization_newton(mesh[ielem], delta_t);
      }

      // Добавляем в глобальную матрицу
      global_A->add_to_elem(0, ielem, local_A[0][0]);
      global_A->add_to_elem(0, ielem + 1, local_A[1][1]);
      global_A->add_to_elem(1, ielem, local_A[1][0]);
      global_A->add_to_elem(2, ielem, local_A[0][1]);

      // Добавляем в глобальный вектор
      global_f[ielem] += local_f[0];
      global_f[ielem + 1] += local_f[1];
   }
}

// Задание начального условия
void mfe::initial_condition(mesh &mesh)
{
   double t0 = mesh.get_time(0);
   q_prev[0] = u(mesh[0].x, t0);
   for (uint32_t i = 0; i < mesh.get_fe_count(); i++)
      q_prev[i + 1] = u(mesh[i].x_next, t0);
}

// Учет первого краевого условия
void mfe::first_boundary_condition(mesh &mesh, const double t, const double delta_t)
{
   // На диагонали ставим 1
   global_A->set_elem(0, 0, 1.0);
   global_A->set_elem(2, 0, 0.0);
   global_A->set_elem(0, mesh.get_fe_count(), 1.0);
   global_A->set_elem(1, mesh.get_fe_count() - 1, 0.0);

   // В вектор правой части ставим известные значения неизвестной функции
   global_f[0] = u(mesh[0].x, t);
   global_f[mesh.get_fe_count()] = u(mesh[mesh.get_fe_count() - 1].x_next, t);

   if (method_to_solve == method::NEWTON)
   {
      A->set_elem(0, 0, 1.0);
      A->set_elem(2, 0, 0.0);
      A->set_elem(0, mesh.get_fe_count(), 1.0);
      A->set_elem(1, mesh.get_fe_count() - 1, 0.0);

      nl_f[0] = u(mesh[0].x, t);
      nl_f[mesh.get_fe_count()] = u(mesh[mesh.get_fe_count() - 1].x_next, t);
   }
}

// Решить систему нелинейных уравнений
std::pair<uint32_t, double> mfe::solve(mesh &mesh, uint32_t max_iter, double eps)
{
   if (method_to_solve == method::NEWTON)
   {
      A = new matrix(mesh.get_fe_count() + 1, 3);
      nl_f.resize(mesh.get_fe_count() + 1);
   }

   initial_condition(mesh);

   solver slv;

   uint32_t iter = 0;

   // Решаем на каждом временном слое
   for (uint32_t i = 1; i < mesh.get_time_size(); i++)
   {
      double delta_t = mesh.get_time(i) - mesh.get_time(i - 1);

      q = q_prev;

      assembly_global_slae(mesh, mesh.get_time(i), delta_t);
      first_boundary_condition(mesh, mesh.get_time(i), delta_t);

      for (; iter < max_iter; iter++)
      {
         q_prev_iter = q;

         slv.solve_by_LU(*global_A, global_f, q);

         assembly_global_slae(mesh, mesh.get_time(i), delta_t);
         first_boundary_condition(mesh, mesh.get_time(i), delta_t);

         double res = 0.0;

         if (residual() < eps)
            break;

         if (use_relax && step() < eps / 2.0)
         {
            double w = 0.0;

            if (method_to_solve == method::SIMPLE_ITERATION)
               w = one_dimensional_search::minimize(*global_A, q, q_prev_iter, global_f, eps);
            else
               w = one_dimensional_search::minimize(*A, q, q_prev_iter, nl_f, eps);

            add_relax(w);

            if (residual() < eps)
               break;
         }
      }
      q_prev = q;
   }

   auto res = std::make_pair(iter, error(mesh));

   //save(directory, res);

   return res;
}

// Добавка параметра релаксации
void mfe::add_relax(const double w)
{
   for (uint32_t i = 0; i < q.size(); i++)
      q[i] = w * q[i] + (1 - w) * q_prev_iter[i];
}

// Разница между решением на прошлой и на текущей итерации
double mfe::step()
{
   std::vector<double> tmp = q;
   for (uint32_t i = 0; i < q.size(); i++)
      tmp[i] -= q_prev_iter[i];

   return norm(tmp) / norm(q);
}

// Расчет невязки
double mfe::residual()
{
   std::unique_ptr<std::vector<double>> Aq;
   if (method_to_solve == method::SIMPLE_ITERATION)
   {
      Aq = global_A->dot(q);
      for (uint32_t i = 0; i < q.size(); i++)
         (*Aq)[i] -= global_f[i];

      return norm(*Aq) / norm(global_f);
   }

   Aq = A->dot(q);
   for (uint32_t i = 0; i < q.size(); i++)
      (*Aq)[i] -= nl_f[i];

   return norm(*Aq) / norm(nl_f);
}

// Расчет относительной невязки
double mfe::error(mesh &mesh)
{
   double err = 0.0;
   double exact_norm = 0.0;
   double t_final = mesh.get_time(mesh.get_time_size() - 1);

   exact.clear();
   exact.resize(q.size(), 0);

   for (uint32_t ielem = 0; ielem < mesh.get_fe_count(); ielem++)
   {
      double h = mesh[ielem].x_next - mesh[ielem].x;
      double xl = mesh[ielem].x;
      double xr = mesh[ielem].x_next;

      // 2-точечная квадратура Гаусса
      double gauss_points[2] = { -0.5773502691896257, 0.5773502691896257 };
      double weights[2] = { 1.0, 1.0 };

      for (int gp = 0; gp < 2; gp++)
      {
         double xi = (gauss_points[gp] + 1.0) / 2.0;
         double x = xl + xi * h;

         // Линейная интерполяция
         double u_num = (1.0 - xi) * q[ielem] + xi * q[ielem + 1];
         double u_exact = u(x, t_final);

         double weight = weights[gp] * h / 2.0;
         err += weight * (u_exact - u_num) * (u_exact - u_num);
         exact_norm += weight * u_exact * u_exact;
      }
   }

   return sqrt(err) / sqrt(exact_norm);
}

// Сохранить результат в файл
void mfe::save(std::string dir, std::pair<uint32_t, double> &result)
{
   std::ofstream res(dir + "result.txt");

   if (res.is_open())
   {
      res << std::left
         << std::setw(15) << "Iterations: " << result.first << std::endl
         << std::setw(15) << "Error: " << result.second << std::endl << std::endl;

      res << std::left << std::setw(15) << "q*"
         << std::setw(15) << "q"
         << std::setw(15) << "|q*-q|" << std::endl;
      res << "------------------------------------------" << std::endl;

      for (uint32_t i = 0; i < q.size(); i++)
      {
         res << std::left << std::setw(15) << exact[i]
            << std::setw(15) << q[i]
            << std::setw(15) << abs(exact[i] - q[i]) << std::endl;
      }

      exact.clear();
      res.close();
   }
}

//---------------------------------------------------------------------------------
//----------------------------- РЕАЛИЗАЦИЯ МЕТОДА НЬЮТОНА -------------------------
//
// Добавки к локальным матрицам и локальному вектору правой части
void mfe::linearization_newton(const finite_elem &elem, const double delta_t)
{
   double h = elem.x_next - elem.x;
   uint32_t ielem = elem.number;

   double dlambda_left = dlambda(q[ielem + 1], q[ielem], h, elem.x, 1);
   double dlambda_right = dlambda(q[ielem + 1], q[ielem], h, elem.x_next, 2);

   double dA_dq1_left = -dlambda_left / h;
   double dA_dq1_right = -dlambda_right / h;
   double dA_dq2_left = dlambda_left / h;
   double dA_dq2_right = dlambda_right / h;

   double dA11_dq1, dA11_dq2, dA12_dq1, dA12_dq2, dA21_dq1, dA21_dq2, dA22_dq1, dA22_dq2;

   dA11_dq1 = dA22_dq1 = (dA_dq1_left + dA_dq1_right) / (2.0 * h);
   dA11_dq2 = dA22_dq2 = (dA_dq2_left + dA_dq2_right) / (2.0 * h);
   dA12_dq1 = dA21_dq1 = -(dA_dq1_left + dA_dq1_right) / (2.0 * h);
   dA12_dq2 = dA21_dq2 = -(dA_dq2_left + dA_dq2_right) / (2.0 * h);

   local_A[0][0] += dA11_dq1 * q[ielem] + dA12_dq1 * q[ielem + 1];
   local_A[0][1] += dA11_dq2 * q[ielem] + dA12_dq2 * q[ielem + 1];
   local_A[1][0] += dA21_dq1 * q[ielem] + dA22_dq1 * q[ielem + 1];
   local_A[1][1] += dA21_dq2 * q[ielem] + dA22_dq2 * q[ielem + 1];

   local_f[0] += q[ielem] * (dA11_dq1 * q[ielem] + dA11_dq2 * q[ielem + 1]) +
      q[ielem + 1] * (dA12_dq1 * q[ielem] + dA12_dq2 * q[ielem + 1]);
   local_f[1] += q[ielem] * (dA21_dq1 * q[ielem] + dA21_dq2 * q[ielem + 1]) +
      q[ielem + 1] * (dA22_dq1 * q[ielem] + dA22_dq2 * q[ielem + 1]);
}

// Производная от функции lambda по qi
double mfe::dlambda(double q2, double q1, double h, double x, uint32_t var)
{
   double du_dx = (q2 - q1) / h;
   //return 0;
   return 1.0;
   //return 2 * du_dx;
   //return exp(du_dx);
}
