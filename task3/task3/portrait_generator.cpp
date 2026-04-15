// ================  PORTRAIT_GENERATOR.CPP ================
#include "portrait_generator.hpp"

/*
  Строит портрет разреженной матрицы размера 2N×2N (N = число узлов).
  Каждый узел порождает пару строк (sin- и cos-компоненты гармонической задачи).

  Структура ig/jg:
    Строка 2i  (sin): связи 2j, 2j+1 для всех j ∈ list[i]  (нижний треугольник)
    Строка 2i+1(cos): связи 2j, 2j+1 для всех j ∈ list[i], плюс 2i (блочная связь)
*/
void portrait_generator::portrait(space_grid &grid,
   std::vector<uint32_t> &ig,
   std::vector<uint32_t> &jg)
{
   // Для каждого узла собираем множество узлов с меньшим номером,
   // с которыми он связан через общий элемент
   std::vector<std::set<uint32_t>> list(grid.get_nodes_count());

   for (uint32_t ielem = 0; ielem < grid.get_elems_count(); ielem++)
   {
      finite_elem elem = grid.get_elem(ielem);

      // 4-узловой билинейный элемент: перебираем все пары (i, j), i < j
      for (uint32_t i = 0; i < 3; i++)
      {
         uint32_t ni = elem[i];
         for (uint32_t j = i + 1; j < 4; j++)
         {
            uint32_t nj = elem[j];
            // nj > ni всегда при нашей нумерации, добавляем ni в list[nj]
            list[nj].insert(ni);
         }
      }
   }

   // Формируем ig (размер 2N+1)
   ig.resize(2 * grid.get_nodes_count() + 1);
   ig[0] = 0;
   ig[1] = 0; // строка 0 (sin) — нет записей ниже диагонали
   ig[2] = 1; // строка 1 (cos) — одна запись: связь с sin того же узла

   for (uint32_t i = 1; i < list.size(); i++)
   {
      uint32_t cnt = static_cast<uint32_t>(list[i].size());
      ig[2 * i + 1] = ig[2 * i] + cnt * 2;       // строка 2i  (sin)
      ig[2 * (i + 1)] = ig[2 * i + 1] + cnt * 2 + 1;   // строка 2i+1(cos)
   }

   // Формируем jg
   jg.resize(ig.back());

   for (uint32_t i = 1, k = 1; i < list.size(); i++)
   {
      // Строка 2i (sin-компонента): столбцы 2j, 2j+1 для j ∈ list[i]
      for (auto j : list[i])
      {
         jg[k] = 2 * j;
         jg[k + 1] = 2 * j + 1;
         k += 2;
      }

      // Строка 2i+1 (cos-компонента): те же столбцы + столбец 2i
      for (auto j : list[i])
      {
         jg[k] = 2 * j;
         jg[k + 1] = 2 * j + 1;
         k += 2;
      }
      jg[k] = 2 * i; // связь cos-строки с sin-строкой того же узла
      k++;
   }
}