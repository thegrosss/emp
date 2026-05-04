// ================  MESH.HPP ================
#include "mesh.hpp"
#include "json.hpp"

// =========================================================
void space_grid::save(std::string path)
{
   std::ofstream grid_out(path + "mesh.json");
   nlohmann::json grid{};

   if (grid_out.is_open())
   {
      for (uint32_t i = 0; i < elems.size(); i++)
         grid["elems"].push_back(elems[i].nodes);

      for (uint32_t i = 0; i < faces.size(); i++)
         grid["faces"].push_back(faces[i].nodes);

      for (uint32_t i = 0; i < points.size(); i++)
         grid["points"].push_back({ points[i].x, points[i].y, points[i].z });

      grid_out << grid << std::endl;

      grid_out.close();
   }
   else throw "Can't open file\n";
}