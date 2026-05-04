// ================  MESH_GENERATOR.HPP ================
#include "mesh_generator.hpp"
#include "json.hpp"

// ===================================================================
space_grid_generator::space_grid_generator()
{
   area_xyz = new area();
   nx = ny = nz = 0;
   kx = ky = kz = 0.0;
   type = mesh_type::UNIFORM;
   nested = 0;
}

void space_grid_generator::read_data(std::string path)
{
   std::ifstream input(path + "space_grid.json");

   if (input.is_open())
   {
      nlohmann::json space_grid{};

      input >> space_grid;

      input.close();

      area_xyz->x_start = space_grid["area"]["x_start"];
      area_xyz->x_end = space_grid["area"]["x_end"];
      area_xyz->y_start = space_grid["area"]["y_start"];
      area_xyz->y_end = space_grid["area"]["y_end"];
      area_xyz->z_start = space_grid["area"]["z_start"];
      area_xyz->z_end = space_grid["area"]["z_end"];
      area_xyz->lambda = space_grid["area"]["lambda"];
      area_xyz->sigma = space_grid["area"]["sigma"];
      area_xyz->hi = space_grid["area"]["hi"];

      nx = space_grid["parameters"]["nx"];
      ny = space_grid["parameters"]["ny"];
      nz = space_grid["parameters"]["nz"];
      kx = space_grid["parameters"]["kx"];
      ky = space_grid["parameters"]["ky"];
      kz = space_grid["parameters"]["kz"];
      type = space_grid["parameters"]["type"];
      nested = space_grid["parameters"]["nested"];
   }
   else throw "Can't open file\n";
}

void space_grid_generator::generate_nodes()
{
   if (nested == 1)
   {
      nx *= 2; kx = sqrt(kx);
      ny *= 2; ky = sqrt(ky);
      nz *= 2; kz = sqrt(kz);
   }
   else if (nested == 2)
   {
      nx *= 4; kx = sqrt(sqrt(kx));
      ny *= 4; ky = sqrt(sqrt(ky));
      nz *= 4; kz = sqrt(sqrt(kz));
   }
   else if (nested == 3)
   {
      nx *= 8; kx = sqrt(sqrt(sqrt(kx)));
      ny *= 8; ky = sqrt(sqrt(sqrt(ky)));
      nz *= 8; kz = sqrt(sqrt(sqrt(kz)));
   }


   if (type == mesh_type::UNIFORM)
   {
      double hx = (area_xyz->x_end - area_xyz->x_start) / double(nx);
      double hy = (area_xyz->y_end - area_xyz->y_start) / double(ny);
      double hz = (area_xyz->z_end - area_xyz->z_start) / double(nz);

      x.resize(nx + 1);
      y.resize(ny + 1);
      z.resize(nz + 1);

      for (uint32_t i = 0; i < x.size(); i++)
         x[i] = area_xyz->x_start + hx * i;

      for (uint32_t i = 0; i < y.size(); i++)
         y[i] = area_xyz->y_start + hy * i;

      for (uint32_t i = 0; i < z.size(); i++)
         z[i] = area_xyz->z_start + hz * i;
   }
   else
   {
      double hx = (area_xyz->x_end - area_xyz->x_start) * (1 - kx) / (1 - pow(kx, nx));
      double hy = (area_xyz->y_end - area_xyz->y_start) * (1 - ky) / (1 - pow(ky, ny));
      double hz = (area_xyz->z_end - area_xyz->z_start) * (1 - kz) / (1 - pow(kz, nz));

      x.resize(nx + 1);
      y.resize(ny + 1);
      z.resize(nz + 1);

      for (uint32_t i = 0; i < x.size(); i++)
         x[i] = area_xyz->x_start + hx * kx * i;

      for (uint32_t i = 0; i < y.size(); i++)
         y[i] = area_xyz->y_start + hy * ky * i;

      for (uint32_t i = 0; i < z.size(); i++)
         z[i] = area_xyz->z_start + hz * kz * i;
   }
}

void space_grid_generator::make_bc(space_grid *&grid, function3D &us, function3D &uc)
{
   // Формируем 1-ые краевые условия ------------------------------
   std::set<uint32_t> dirichlet;

   uint32_t dirichlet_counter = 0;

   for (uint32_t i = 0; i < grid->faces.size(); i++)
      for (uint32_t j = 0; j < grid->faces[i].nodes.size(); j++)
         dirichlet.insert(grid->faces[i].nodes[j]);

   grid->dirichlet.resize(2 * dirichlet.size());

   uint32_t i = 0;

   for (const auto &it : dirichlet)
   {
      point3D point = grid->get_point(it);
      grid->dirichlet[2 * i] = { 2 * it, us(point.x, point.y, point.z) };
      grid->dirichlet[2 * i + 1] = { 2 * it + 1, uc(point.x, point.y, point.z) };
      i++;
   }

   dirichlet.clear();
}

void space_grid_generator::build_mesh(space_grid *&grid, function3D &us, function3D &uc)
{
   read_data();
   generate_nodes();

   uint32_t n_xy_points = (nx + 1) * (ny + 1);

   grid = new space_grid(nx, ny, nz);

   grid->set_type(type);

   // Формируем узлы -------------------------------------------
   uint32_t pos = 0;

   for (uint32_t i = 0; i < z.size(); i++)
      for (uint32_t j = 0; j < y.size(); j++)
         for (uint32_t k = 0; k < x.size(); k++)
            grid->points[pos++] = point3D(x[k], y[j], z[i]);

   x.clear();
   y.clear();
   z.clear();

   // Формируем конечные элементы ------------------------------
   std::array<uint32_t, 8> nodes;

   uint32_t ielem = 0;

   for (uint32_t k = 0; k < nz; k++)
   {
      for (uint32_t i = 0; i < ny; i++)
      {
         for (uint32_t j = 0; j < nx; j++)
         {
            nodes[0] = j + i * (nx + 1) + k * n_xy_points;
            nodes[1] = j + i * (nx + 1) + 1 + k * n_xy_points;
            nodes[2] = j + i * (nx + 1) + nx + 1 + k * n_xy_points;
            nodes[3] = j + i * (nx + 1) + nx + 2 + k * n_xy_points;
            nodes[4] = j + i * (nx + 1) + n_xy_points + k * n_xy_points;
            nodes[5] = j + i * (nx + 1) + n_xy_points + 1 + k * n_xy_points;
            nodes[6] = j + i * (nx + 1) + n_xy_points + nx + 1 + k * n_xy_points;
            nodes[7] = j + i * (nx + 1) + n_xy_points + nx + 2 + k * n_xy_points;

            grid->elems[ielem++] = finite_elem(
               nodes, area_xyz->lambda, area_xyz->sigma, area_xyz->hi
            );
         }
      }
   }

   delete area_xyz;
   area_xyz = nullptr;

   // Формируем грани ------------------------------------------
   std::array<uint32_t, 4> face_nodes;

   // XY при Z = Z_min и Z = Z_max
   uint32_t iface = 0;

   for (uint32_t k = 0; k < 2; k++)
   {
      for (uint32_t i = 0; i < ny; i++)
      {
         for (uint32_t j = 0; j < nx; j++)
         {
            face_nodes[0] = n_xy_points * nz * k + j + i * (nx + 1);
            face_nodes[1] = n_xy_points * nz * k + j + i * (nx + 1) + 1;
            face_nodes[2] = n_xy_points * nz * k + j + i * (nx + 1) + nx + 1;
            face_nodes[3] = n_xy_points * nz * k + j + i * (nx + 1) + nx + 2;

            grid->faces[iface++] = face(face_nodes, iface);
         }
      }
   }

   // XZ при Y = Y_min и Y = Y_max
   for (uint32_t k = 0; k < 2; k++)
   {
      for (uint32_t i = 0; i < nz; i++)
      {
         for (uint32_t j = 0; j < nx; j++)
         {
            face_nodes[0] = k * (nx + 1) * ny + j + i * n_xy_points;
            face_nodes[1] = k * (nx + 1) * ny + j + i * n_xy_points + 1;
            face_nodes[2] = k * (nx + 1) * ny + j + i * n_xy_points + n_xy_points;
            face_nodes[3] = k * (nx + 1) * ny + j + i * n_xy_points + n_xy_points + 1;

            grid->faces[iface++] = face(face_nodes, iface);
         }
      }
   }

   // YZ при X = X_min и X = X_max
   for (uint32_t k = 0; k < 2; k++)
   {
      for (uint32_t i = 0; i < nz; i++)
      {
         for (uint32_t j = 0; j < ny; j++)
         {
            face_nodes[0] = k * nx + j * (nx + 1) + i * n_xy_points;
            face_nodes[1] = k * nx + j * (nx + 1) + i * n_xy_points + nx + 1;
            face_nodes[2] = k * nx + j * (nx + 1) + i * n_xy_points + n_xy_points;
            face_nodes[3] = k * nx + j * (nx + 1) + i * n_xy_points + n_xy_points + nx + 1;

            grid->faces[iface++] = face(face_nodes, iface);
         }
      }
   }

   make_bc(grid, us, uc);
}