// ================  MFE.CPP ================
#include "mfe.hpp"

// ----------------------------------------------------------------
mfe::mfe(space_grid &mesh)
{
   grid = &mesh;
   w = 0.0;

   A = new sparse_matrix(2 * grid->get_nodes_count());

   // Инициализируем 4×4 локальные матрицы (билинейные базисные функции)
   local_p.assign(4, dvector(4, 0.0));
   local_c.assign(4, dvector(4, 0.0));
   G.assign(4, dvector(4, 0.0));
   M.assign(4, dvector(4, 0.0));

   local_vec_fs.resize(4);
   local_vec_fc.resize(4);
}

// ----------------------------------------------------------------
// Сборка локальных матриц G (жёсткости) и M (масс),
// затем P = λ*G - ω²χ*M  и  C = ω*σ*M
void mfe::build_local_matrices(finite_elem &elem)
{
   // Опорный элемент: node[0] — нижний левый, node[3] — верхний правый
   omega om = { grid->get_point(elem[0]), grid->get_point(elem[3]) };
   basis_function bf(om);

   function2D f;

   for (uint8_t i = 0; i < 4; i++)
   {
      for (uint8_t j = 0; j <= i; j++)
      {
         // G_ij = ∫(∇ψ_i · ∇ψ_j) dΩ
         f = [&](double x, double y)
            {
               point2D pt(x, y);
               return bf.d_psi(i, 1, pt) * bf.d_psi(j, 1, pt)
                  + bf.d_psi(i, 2, pt) * bf.d_psi(j, 2, pt);
            };
         G[i][j] = G[j][i] = gauss.integrate2D(f, om);

         // M_ij = ∫(ψ_i · ψ_j) dΩ
         f = [&](double x, double y)
            {
               point2D pt(x, y);
               return bf.psi(i, pt) * bf.psi(j, pt);
            };
         M[i][j] = M[j][i] = gauss.integrate2D(f, om);

         local_p[i][j] = local_p[j][i] = elem.lambda * G[i][j] - w * w * elem.hi * M[i][j];
         local_c[i][j] = local_c[j][i] = w * elem.sigma * M[i][j];
      }
   }
}

// ----------------------------------------------------------------
// Сборка локальных векторов правой части через матрицу масс M
void mfe::build_local_vectors(finite_elem &elem, r_function2D &fs, r_function2D &fc)
{
   dvector fs_vec(4), fc_vec(4);
   local_vec_fs.assign(4, 0.0);
   local_vec_fc.assign(4, 0.0);

   for (uint8_t i = 0; i < 4; i++)
   {
      point2D pt = grid->get_point(elem[i]);
      fs_vec[i] = fs(pt.x, pt.y, elem.lambda, elem.sigma, elem.hi, w);
      fc_vec[i] = fc(pt.x, pt.y, elem.lambda, elem.sigma, elem.hi, w);
   }

   for (uint8_t i = 0; i < 4; i++)
      for (uint8_t j = 0; j < 4; j++)
      {
         local_vec_fs[i] += M[i][j] * fs_vec[j];
         local_vec_fc[i] += M[i][j] * fc_vec[j];
      }
}

// ----------------------------------------------------------------
void mfe::add_to_global_matrix(const uint32_t i, const uint32_t j, const double val)
{
   if (i == j)
   {
      A->di[i] += val;
      return;
   }

   if (i < j)
   {
      for (uint32_t ind = A->ig[j]; ind < A->ig[j + 1]; ind++)
         if (A->jg[ind] == i) { A->ggu[ind] += val; return; }
   }
   else
   {
      for (uint32_t ind = A->ig[i]; ind < A->ig[i + 1]; ind++)
         if (A->jg[ind] == j) { A->ggl[ind] += val; return; }
   }
}

// ----------------------------------------------------------------
/*
  Глобальная матрица гармонической задачи имеет блочную структуру 2×2:
      | P  -C |  | q_s |   | f_s |
      |       |  |     | = |     |
      | C   P |  | q_c |   | f_c |

  Добавляем вклад каждого элемента: P_ij в (2i,2j) и (2i+1,2j+1),
                                     C_ij в (2i+1,2j) и -C_ij в (2i,2j+1).
*/
void mfe::assembly_global_matrix_and_vector(r_function2D &fs, r_function2D &fc)
{
   portrait_generator::portrait(*grid, A->ig, A->jg);

   A->di.resize(A->dim, 0.0);
   A->ggl.resize(A->ig.back(), 0.0);
   A->ggu.resize(A->ig.back(), 0.0);
   b.resize(A->dim, 0.0);

   for (uint32_t ielem = 0; ielem < grid->get_elems_count(); ielem++)
   {
      finite_elem elem = grid->get_elem(ielem);

      build_local_matrices(elem);
      build_local_vectors(elem, fs, fc);

      for (uint8_t i = 0; i < 4; i++)
      {
         b[2 * elem[i]] += local_vec_fs[i];
         b[2 * elem[i] + 1] += local_vec_fc[i];

         for (uint8_t j = 0; j < 4; j++)
         {
            add_to_global_matrix(2 * elem[i], 2 * elem[j], local_p[i][j]);
            add_to_global_matrix(2 * elem[i] + 1, 2 * elem[j] + 1, local_p[i][j]);
            add_to_global_matrix(2 * elem[i] + 1, 2 * elem[j], local_c[i][j]);
            add_to_global_matrix(2 * elem[i], 2 * elem[j] + 1, -local_c[i][j]);
         }
      }
   }

   add_dirichlet();
}

// ----------------------------------------------------------------
void mfe::add_dirichlet()
{
   for (uint32_t i = 0; i < grid->get_first_bound_count(); i++)
   {
      auto cond = grid->get_dirichlet_cond(i);

      A->di[cond.node] = 1.0;
      b[cond.node] = cond.value;

      // Обнуляем нижнюю строку (ggl)
      for (uint32_t k = A->ig[cond.node]; k < A->ig[cond.node + 1]; k++)
         A->ggl[k] = 0.0;

      // Обнуляем верхнюю строку (ggu) во всех зависимых строках
      for (uint32_t k = cond.node + 1; k < A->size(); k++)
         for (uint32_t j = A->ig[k]; j < A->ig[k + 1]; j++)
            if (A->jg[j] == cond.node) { A->ggu[j] = 0.0; break; }
   }
}

// ----------------------------------------------------------------
result mfe::solve(method method_to_solve)
{
   q.clear();

   solver slv;

   if (method_to_solve == method::LU)
      return slv.solve_by_LU(*(A->to_profile()), b, q);

   // LOS_LU — локально-оптимальная схема с LU-предобусловливанием
   return slv.solve_iterative(*A, b, q, method::LOS_LU, 1000, 1e-16);
}