// ================  PORTRAIT_GENERATOR.HPP ================
#include "portrait_generator.hpp"

void portrait_generator::portrait(space_grid &grid, std::vector<uint32_t> &ig, std::vector<uint32_t> &jg)
{
   std::vector<std::set<uint32_t>> list(grid.get_nodes_count());

   for (uint32_t ielem = 0; ielem < grid.get_elems_count(); ielem++)
   {
      finite_elem elem = grid.get_elem(ielem);

      for (uint32_t i = 0; i < 7; i++)
      {
         uint32_t elem_to_insert = elem[i];

         for (uint32_t j = i + 1; j < 8; j++)
         {
            uint32_t insert_pos = elem[j];

            list[insert_pos].insert(elem_to_insert);
         }
      }
   }

   ig.resize(2 * grid.get_nodes_count() + 1);
   ig[0] = 0;
   ig[1] = 0;
   ig[2] = 1;

   for (uint32_t i = 1; i < list.size(); i++)
   {
      ig[2 * i + 1] = ig[2 * i] + list[i].size() * 2;
      ig[2 * (i + 1)] = ig[2 * i + 1] + list[i].size() * 2 + 1;
   }

   jg.resize(ig.back());

   for (uint32_t i = 1, k = 1; i < list.size(); i++)
   {
      for (auto j : list[i])
      {
         jg[k] = 2 * j;
         jg[k + 1] = 2 * j + 1;
         k += 2;
      }
      for (auto j : list[i])
      {
         jg[k] = 2 * j;
         jg[k + 1] = 2 * j + 1;
         k += 2;
      }
      jg[k] = 2 * i;
      k++;
   }
}