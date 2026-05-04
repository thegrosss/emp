// ================  MFE.CPP ================
#include "mfe.hpp"

mfe::mfe(space_grid &mesh)
{
   grid = &mesh;

   w = 0.0;

   A = new sparse_matrix(2 * grid->get_nodes_count());

   local_p.resize(8);
   local_c.resize(8);

   G.resize(8);
   M.resize(8);

   local_vec_fs.resize(8);
   local_vec_fc.resize(8);

   for (uint8_t i = 0; i < 8; i++)
   {
      local_p[i].resize(8);
      local_c[i].resize(8);

      G[i].resize(8);
      M[i].resize(8);
   }
}

void mfe::build_local_matrices(finite_elem &elem)
{
   omega omega = { point3D(grid->get_point(elem[0])), point3D(grid->get_point(elem[7])) };

   basis_function bf(omega);

   function3D f;

   for (uint8_t i = 0; i < 8; i++)
   {
      for (uint8_t j = 0; j <= i; j++)
      {
         f = [&](double x, double y, double z)
            {
               point3D point(x, y, z);

               double d_psi_i_x = bf.d_psi(i, 1, point);
               double d_psi_j_x = bf.d_psi(j, 1, point);

               double d_psi_i_y = bf.d_psi(i, 2, point);
               double d_psi_j_y = bf.d_psi(j, 2, point);

               double d_psi_i_z = bf.d_psi(i, 3, point);
               double d_psi_j_z = bf.d_psi(j, 3, point);

               return (d_psi_i_x * d_psi_j_x + d_psi_i_y * d_psi_j_y + d_psi_i_z * d_psi_j_z);
            };

         G[i][j] = G[j][i] = gauss.integrate3D(f, omega);


         // -----------------------------------------------------------
         f = [&](double x, double y, double z)
            {
               point3D point(x, y, z);

               double psi_i = bf.psi(i, point);
               double psi_j = bf.psi(j, point);

               return psi_i * psi_j;
            };

         M[i][j] = M[j][i] = gauss.integrate3D(f, omega);

         local_p[i][j] = local_p[j][i] = elem.lambda * G[i][j] - w * w * elem.hi * M[i][j];
         local_c[i][j] = local_c[j][i] = w * elem.sigma * M[i][j];
      }
   }
}

void mfe::build_local_vectors(finite_elem &elem, r_function3D &fs, r_function3D &fc)
{
   dvector fs_vec(8);
   dvector fc_vec(8);

   local_vec_fs.clear();
   local_vec_fc.clear();

   local_vec_fs.resize(8);
   local_vec_fc.resize(8);


   for (uint8_t i = 0; i < 8; i++)
   {
      point3D point = grid->get_point(elem[i]);

      fs_vec[i] = fs(point.x, point.y, point.z, elem.lambda, elem.sigma, elem.hi, w);
      fc_vec[i] = fc(point.x, point.y, point.z, elem.lambda, elem.sigma, elem.hi, w);
   }

   for (uint8_t i = 0; i < 8; i++)
   {
      for (uint8_t j = 0; j < 8; j++)
      {
         local_vec_fs[i] += M[i][j] * fs_vec[j];
         local_vec_fc[i] += M[i][j] * fc_vec[j];
      }
   }

   fs_vec.clear();
   fc_vec.clear();
}

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
      {
         if (A->jg[ind] == i)
         {
            A->ggu[ind] += val;
            return;
         }
      }
   }
   else
   {
      for (uint32_t ind = A->ig[i]; ind < A->ig[i + 1]; ind++)
      {
         if (A->jg[ind] == j)
         {
            A->ggl[ind] += val;

         }
      }
   }
}

void mfe::assembly_global_matrix_and_vector(r_function3D &fs, r_function3D &fc)
{
   portrait_generator::portrait(*grid, A->ig, A->jg);

   A->di.resize(A->dim);
   A->ggl.resize(A->ig.back());
   A->ggu.resize(A->ig.back());
   b.resize(A->dim);

   for (uint32_t ielem = 0; ielem < grid->get_elems_count(); ielem++)
   {
      finite_elem elem = grid->get_elem(ielem);

      build_local_matrices(elem);

      build_local_vectors(elem, fs, fc);

      for (uint8_t i = 0; i < 8; i++)
      {
         b[2 * elem[i]] += local_vec_fs[i];
         b[2 * elem[i] + 1] += local_vec_fc[i];

         for (uint8_t j = 0; j < 8; j++)
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

void mfe::add_dirichlet()
{
   for (uint32_t i = 0; i < grid->get_first_bound_count(); i++)
   {
      auto cond = grid->get_dirichlet_cond(i);

      A->di[cond.node] = 1.0;
      b[cond.node] = cond.value;

      // Обнуляем строку
      uint32_t i0 = A->ig[cond.node];
      uint32_t i1 = A->ig[cond.node + 1];

      for (uint32_t k = i0; k < i1; k++)
         A->ggl[k] = 0.0;

      for (uint32_t k = i + 1; k < A->size(); k++)
      {
         uint32_t i0 = A->ig[k];
         uint32_t i1 = A->ig[k + 1];

         for (uint32_t j = i0; j < i1; j++)
         {
            if (A->jg[j] == cond.node)
               A->ggu[j] = 0.0;
         }
      }
   }
}

result mfe::solve(method method_to_solve)
{
   q.clear();

   result res;

   res.iters = 0;
   res.residual = 0;
   res.time = std::chrono::microseconds(0);

   solver slv;

   if (method_to_solve == method::LU)
   {
      res = slv.solve_by_LU(*(A->to_profile()), b, q);
   }
   else if (method_to_solve == method::LOS_LU)
   {
      res = slv.solve_iterative(*A, b, q, method::LOS_LU, 1000, 1e-16);
   }
   else if (method_to_solve == method::BCGSTAB_LU)
   {
      res = slv.solve_iterative(*A, b, q, method::BCGSTAB_LU, 1000, 1e-16);
   }

   return res;
}