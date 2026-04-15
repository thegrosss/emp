// ================  MESH_GENERATOR.CPP ================
#include "mesh_generator.hpp"

// ====================================================================
space_grid_generator::space_grid_generator()
{
   area_xy = new area();
   nx = ny = 0;
   kx = ky = 0.0;
   type = mesh_type::UNIFORM;
   nested = 0;
}

// ----------------------------------------------------------------
void space_grid_generator::read_data(std::string path)
{
   std::ifstream input(path + "space_grid.json");

   if (input.is_open())
   {
      nlohmann::json sg{};
      input >> sg;
      input.close();

      area_xy->x_start = sg["area"]["x_start"];
      area_xy->x_end = sg["area"]["x_end"];
      area_xy->y_start = sg["area"]["y_start"];
      area_xy->y_end = sg["area"]["y_end"];
      area_xy->lambda = sg["area"]["lambda"];
      area_xy->sigma = sg["area"]["sigma"];
      area_xy->hi = sg["area"]["hi"];

      nx = sg["parameters"]["nx"];
      ny = sg["parameters"]["ny"];
      kx = sg["parameters"]["kx"];
      ky = sg["parameters"]["ky"];
      type = sg["parameters"]["type"];
      nested = sg["parameters"]["nested"];
   }
   else throw "Can't open space_grid.json\n";
}

// ----------------------------------------------------------------
void space_grid_generator::generate_nodes()
{
   // Вложенные измельчения сетки
   if (nested == 1) { nx *= 2; kx = sqrt(kx); ny *= 2; ky = sqrt(ky); }
   else if (nested == 2) { nx *= 4; kx = sqrt(sqrt(kx)); ny *= 4; ky = sqrt(sqrt(ky)); }
   else if (nested == 3) { nx *= 8; kx = sqrt(sqrt(sqrt(kx))); ny *= 8; ky = sqrt(sqrt(sqrt(ky))); }

   x.resize(nx + 1);
   y.resize(ny + 1);

   if (type == mesh_type::UNIFORM)
   {
      double hx = (area_xy->x_end - area_xy->x_start) / double(nx);
      double hy = (area_xy->y_end - area_xy->y_start) / double(ny);

      for (uint32_t i = 0; i <= nx; i++) x[i] = area_xy->x_start + hx * i;
      for (uint32_t i = 0; i <= ny; i++) y[i] = area_xy->y_start + hy * i;
   }
   else // неравномерная сетка (геометрическая прогрессия)
   {
      double hx = (area_xy->x_end - area_xy->x_start) * (1.0 - kx) / (1.0 - pow(kx, nx));
      double hy = (area_xy->y_end - area_xy->y_start) * (1.0 - ky) / (1.0 - pow(ky, ny));

      x[0] = area_xy->x_start;
      for (uint32_t i = 1; i <= nx; i++) x[i] = x[i - 1] + hx * pow(kx, i - 1);

      y[0] = area_xy->y_start;
      for (uint32_t i = 1; i <= ny; i++) y[i] = y[i - 1] + hy * pow(ky, i - 1);
   }
}

// ----------------------------------------------------------------
// Наложение первых краевых условий (все граничные узлы — Дирихле)
void space_grid_generator::make_bc(space_grid *&grid, function2D &us, function2D &uc)
{
   std::set<uint32_t> dir_nodes;

   for (uint32_t i = 0; i < grid->edges.size(); i++)
      for (uint32_t j = 0; j < grid->edges[i].nodes.size(); j++)
         dir_nodes.insert(grid->edges[i].nodes[j]);

   grid->dirichlet.resize(2 * dir_nodes.size());

   uint32_t idx = 0;
   for (const auto &n : dir_nodes)
   {
      point2D pt = grid->get_point(n);
      grid->dirichlet[2 * idx] = { 2 * n,     us(pt.x, pt.y) };
      grid->dirichlet[2 * idx + 1] = { 2 * n + 1, uc(pt.x, pt.y) };
      idx++;
   }
}

// ----------------------------------------------------------------
void space_grid_generator::build_mesh(space_grid *&grid, function2D &us, function2D &uc)
{
   read_data();
   generate_nodes();

   grid = new space_grid(nx, ny);
   grid->set_type(type);

   // ----- Заполнение узлов (обход: y — внешний, x — внутренний) -----
   uint32_t pos = 0;
   for (uint32_t i = 0; i <= ny; i++)
      for (uint32_t j = 0; j <= nx; j++)
         grid->points[pos++] = point2D(x[j], y[i]);

   x.clear();
   y.clear();

   // ----- Заполнение конечных элементов ------------------------------
   /*
     Нумерация узлов внутри элемента:
       2 --- 3
       |     |
       0 --- 1
     nodes[k] = глобальный номер k-го узла элемента
   */
   std::array<uint32_t, 4> nodes;
   uint32_t ielem = 0;

   for (uint32_t i = 0; i < ny; i++)
   {
      for (uint32_t j = 0; j < nx; j++)
      {
         nodes[0] = i * (nx + 1) + j;
         nodes[1] = i * (nx + 1) + j + 1;
         nodes[2] = (i + 1) * (nx + 1) + j;
         nodes[3] = (i + 1) * (nx + 1) + j + 1;

         grid->elems[ielem++] = finite_elem(
            nodes, area_xy->lambda, area_xy->sigma, area_xy->hi
         );
      }
   }

   delete area_xy;
   area_xy = nullptr;

   // ----- Граничные рёбра (все 4 стороны прямоугольника) ------------
   std::array<uint32_t, 2> edge_nodes;
   uint32_t iedge = 0;

   // Нижняя сторона (y = y_start, строка i = 0)
   for (uint32_t j = 0; j < nx; j++)
   {
      edge_nodes = { j, j + 1 };
      grid->edges[iedge++] = edge(edge_nodes, iedge);
   }

   // Верхняя сторона (y = y_end, строка i = ny)
   for (uint32_t j = 0; j < nx; j++)
   {
      edge_nodes = { ny * (nx + 1) + j, ny * (nx + 1) + j + 1 };
      grid->edges[iedge++] = edge(edge_nodes, iedge);
   }

   // Левая сторона (x = x_start, столбец j = 0)
   for (uint32_t i = 0; i < ny; i++)
   {
      edge_nodes = { i * (nx + 1), (i + 1) * (nx + 1) };
      grid->edges[iedge++] = edge(edge_nodes, iedge);
   }

   // Правая сторона (x = x_end, столбец j = nx)
   for (uint32_t i = 0; i < ny; i++)
   {
      edge_nodes = { i * (nx + 1) + nx, (i + 1) * (nx + 1) + nx };
      grid->edges[iedge++] = edge(edge_nodes, iedge);
   }

   make_bc(grid, us, uc);
}