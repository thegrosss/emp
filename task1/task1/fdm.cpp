#include "fdm.h"

const double RELAXATION_PARAMETER = 1.49;

// Преобразовываем сетку в СЛАУ
void fdm::mesh_to_slae(mesh &mesh, func2D_u &u, func2D_f &f)
{
   A = new diagonal_matrix(mesh.size(), 5, mesh.get_width() - 1);
   b.resize(mesh.size());

   // Количество элементов по оси X
   uint32_t kx = mesh.get_width() + 1;

   for (uint32_t i = 0; i < mesh.size(); i++)
   {
      auto node = mesh[i];

      // Внутренний узел
      if (node.type == 0)
      {
         // Если сетка равномерная
         if (mesh.get_type() == mesh::mesh_type::UNIFORM)
         {
            double hx = mesh[i + 1].p.x - node.p.x;

            double hy = mesh[i + kx].p.y - node.p.y;

            A->diags[0][i] = (2.0 / (hx * hx) + 2.0 / (hy * hy)) * node.lambda + node.gamma;
            A->diags[1][i - 1] = -node.lambda / (hx * hx);
            A->diags[3][i] = -node.lambda / (hx * hx);
            A->diags[2][i - kx] = -node.lambda / (hy * hy);
            A->diags[4][i] = -node.lambda / (hy * hy);
         }
         // Если сетка неравномерная
         //else
         //{
         //   double hx = mesh[i + 1].p.x - node.p.x;
         //   double hx_prev = node.p.x - mesh[i - 1].p.x;

         //   double hy = mesh[i + kx].p.y - node.p.y;
         //   double hy_prev = node.p.y - mesh[i - kx].p.y;

         //   A->diags[0][i] = (2.0 / (hx_prev * hx) + 2.0 / (hy_prev * hy)) * node.lambda + node.gamma;
         //   A->diags[1][i - 1] = -2.0 * node.lambda / (hx_prev * (hx + hx_prev));
         //   A->diags[3][i] = -2.0 * node.lambda / (hx * (hx + hx_prev));
         //   A->diags[2][i - kx] = -2.0 * node.lambda / (hy_prev * (hy + hy_prev));
         //   A->diags[4][i] = -2.0 * node.lambda / (hy * (hy + hy_prev));
         //}

         //b[i] = f(node.p.x, node.p.y, node.lambda, node.gamma);
      }

      // Фиктивный узел
      else if (node.type == -1)
      {
         A->diags[0][i] = 1.0;
         b[i] = 0.0;
      }

      // Граничный узел
      else if (node.type == 1 || node.type == 2 || node.type == 3 || node.type == 4)
      {
         switch (mesh[i].bc)
         {
            // Первые краевые
         case border::bound_cond::DIRICHLET:
         {
            A->diags[0][i] = 1.0;
            b[i] = u(node.p.x, node.p.y);
            break;
         }

         // Вторые краевые
         case border::bound_cond::NEUMANN:
         {
            double hx, hy;
            if (node.type == 2)
            {
               hx = mesh[i + 1].p.x - node.p.x;
               A->diags[0][i] = 1.0 / hx * node.lambda;
               A->diags[3][i] = -1.0 / hx * node.lambda;
               b[i] = -node.lambda * du_dx(u, node.p, hx);
            }
            else if (node.type == 4)
            {
               hx = node.p.x - mesh[i - 1].p.x;
               A->diags[0][i] = 1.0 / hx * node.lambda;
               A->diags[1][i - 1] = -1.0 / hx * node.lambda;
               b[i] = node.lambda * du_dx(u, node.p, hx);
            }
            else if (node.type == 1)
            {
               hy = mesh[i + kx].p.y - node.p.y;
               A->diags[0][i] = 1.0 / hy * node.lambda;
               A->diags[4][i] = -1.0 / hy * node.lambda;
               b[i] = -node.lambda * du_dy(u, node.p, hy);
            }
            else if (node.type == 3)
            {
               hy = node.p.y - mesh[i - kx].p.y;
               A->diags[0][i] = 1.0 / hy * node.lambda;
               A->diags[2][i - kx] = -1.0 / hy * node.lambda;
               b[i] = node.lambda * du_dy(u, node.p, hy);
            }
            break;
         }

         // По-умолчанию задаем первые краевые
         default:
         {
            A->diags[0][i] = 1.0;
            b[i] = u(mesh[i].p.x, mesh[i].p.y);
            break;
         }
         }
      }
   }
}


// Решить систему и вернуть результат
std::pair<uint32_t, double> fdm::calculate(std::vector<double> &u)
{
   u.resize(b.size());

   solver slv(100000, 1e-13);

   auto res = slv.solve(RELAXATION_PARAMETER, *A, b, u);

   res.second = residual();

   return res;
}

// 1-я производная по x
double fdm::du_dx(func2D_u &f, point &p, double hx)
{
   return (f(p.x + hx, p.y) - f(p.x, p.y)) / hx;
}

// 1-я производная по y
double fdm::du_dy(func2D_u &f, point &p, double hy)
{
   return (f(p.x, p.y + hy) - f(p.x, p.y)) / hy;
}

// Посчитать невязку
double fdm::residual()
{
   double sum = 0.0;
   double norm_b = 0.0;

   auto res = A->dot(exact);

   for (uint32_t i = 0; i < b.size(); i++)
   {
      sum += (b[i] - (*res)[i]) * (b[i] - (*res)[i]);
      norm_b += b[i] * b[i];
   }

   return sqrt(sum) / sqrt(norm_b);
}