#include "matrix.h"

matrix::matrix(const uint32_t _size, const uint32_t _tape_width)
{
   size = _size;
   tape_width = _tape_width;

   uint32_t k = (tape_width - 1) / 2;

   tape.resize(tape_width);

   tape[0].resize(size);
   for (uint32_t i = 1; i < k + 1; i++)
   {
      tape[i].resize(size - i);
      tape[k + i].resize(size - i);
   }
}

// Умножить матрицу на вектор
std::unique_ptr<std::vector<double>> matrix::dot(const std::vector<double> &vector)
{
   auto res = std::make_unique<std::vector<double>>(vector.size());

   uint32_t k = (tape_width - 1) / 2;

   for (uint32_t i = 0; i < k + 1; i++)
   {
      uint32_t offset = i;
      for (uint32_t j = 0; j < tape[i].size(); j++, offset++)
      {
         (*res)[offset] += tape[i][j] * vector[j];

         if (i == 0) continue;

         (*res)[j] += tape[k + i][j] * vector[offset];
      }
   }

   return res;
}

// Перевести матрицу в плотный формат
void matrix::to_dense(const std::string dir)
{
   uint32_t k = (tape_width - 1) / 2;

   mat.resize(size);
   for (uint32_t i = 0; i < size; i++)
      mat[i].resize(size);

   for (uint32_t i = 0; i < k + 1; i++)
   {
      uint32_t offset = i;
      for (uint32_t j = 0; j < tape[i].size(); j++, offset++)
      {
         mat[offset][j] = tape[i][j];

         if (i == 0)
            continue;

         mat[j][offset] = tape[k + i][j];
      }
   }

   //===============================================================
   std::ofstream dense(dir + "dense.txt");

   if (dense.is_open())
   {
      dense.precision(6);
      for (uint32_t i = 0; i < size; i++)
      {
         for (uint32_t j = 0; j < size; j++)
         {
            if (mat[i][j] == 0.0)
               dense << std::setw(10) << std::left << ".";
            else
               dense << std::setw(10) << std::left << mat[i][j];
         }
         dense << std::endl;
      }
      dense.close();
   }
   else
      throw "Can't open file";

   //mat.clear();
}

// Обнулить все элементы
void matrix::set_elem_to_null()
{
   for (uint32_t i = 0; i < tape_width; i++)
      for (uint32_t j = 0; j < tape[i].size(); j++)
         tape[i][j] = 0.0;
}