// ================  MATRIX.HPP ================
#include "matrix.hpp"
#include "json.hpp"

std::unique_ptr<dvector> sparse_matrix::dot(const dvector &vector)
{
   auto result = std::make_unique <dvector>(vector.size());

   for (uint32_t i = 0; i < result->size(); i++)
   {
      (*result)[i] = di[i] * vector[i];
      for (uint32_t j = ig[i]; j < ig[i + 1]; j++)
      {
         (*result)[i] += ggl[j] * vector[jg[j]];
         (*result)[jg[j]] += ggu[j] * vector[i];
      }
   }

   return result;
}

void sparse_matrix::to_dense()
{
   dense_matrix.resize(dim);

   for (uint32_t i = 0; i < dense_matrix.size(); i++)
      dense_matrix[i].resize(dim, 0);

   for (uint32_t i = 0; i < dim; i++)
   {
      dense_matrix[i][i] = di[i];
      for (uint32_t j = ig[i]; j < ig[i + 1]; j++)
      {
         dense_matrix[i][jg[j]] = ggl[j];
         dense_matrix[jg[j]][i] = ggu[j];
      }
   }
}

profile_matrix *sparse_matrix::to_profile()
{
   profile_matrix *matrix = new profile_matrix(dim);

   matrix->di = di;

   matrix->ig.resize(dim + 1);
   matrix->ig[0] = 0;

   uint32_t count = 0;

   // Меняем портрет ----------------------------------------------
   for (uint32_t i = 0; i < dim; i++)
   {
      count = 0;

      uint32_t i0 = ig[i];
      uint32_t i1 = ig[i + 1];

      uint32_t prof = i1 - i0;

      if (prof > 0)
      {
         count = i - jg[i0];

         matrix->ig[i + 1] = matrix->ig[i] + count;
      }
      else
      {
         matrix->ig[i + 1] = matrix->ig[i];
      }
   }

   matrix->ggl.resize(matrix->ig.back());
   matrix->ggu.resize(matrix->ig.back());

   // Заполняем верхний и нижний треугольники ---------------------
   for (uint32_t i = 0; i < dim; i++)
   {
      uint32_t i0_p = matrix->ig[i];
      uint32_t i1_p = matrix->ig[i + 1];

      uint32_t j0_p = i - (i1_p - i0_p);

      uint32_t i0_s = ig[i];

      for (uint32_t row_ind = i0_p; row_ind < i1_p; row_ind++, j0_p++)
      {
         if (j0_p == jg[i0_s])
         {
            matrix->ggl[row_ind] = ggl[i0_s];
            matrix->ggu[row_ind] = ggu[i0_s];
            i0_s++;
         }
         else
         {
            matrix->ggl[row_ind] = 0.0;
            matrix->ggu[row_ind] = 0.0;
         }
      }
   }

   return matrix;
}

void sparse_matrix::save(matrix_type type, std::string path)
{
   switch (type)
   {
   case matrix_type::DENSE: {
      std::ofstream dense_out(path + "dense.txt");
      nlohmann::json dense{};

      dense["dense"] = dense_matrix;

      dense_out.precision(3);
      for (size_t i = 0; i < dim; i++)
      {
         for (size_t j = 0; j < dim; j++)
         {
            //if (i >= j)
            dense_out << std::setw(15) << std::left << dense_matrix[i][j];
            //else
            //	dense_out << std::setw(8) << std::left << "*";
         }
         dense_out << std::endl;
      }

      //dense_out << dense;

      dense_out.close();
   }
   break;

   case matrix_type::PROFILE: {
      std::ofstream profile_out(path + "profile.json");
      nlohmann::json profile{};

      profile["ig"] = ig;
      profile["di"] = di;
      profile["ggl"] = ggl;
      profile["ggu"] = ggu;

      profile_out << profile;
      profile_out.close();
   }
   break;

   case matrix_type::SPARSE: {
      std::ofstream sparse_out(path + "sparse.json");
      nlohmann::json sparse{};

      sparse["ig"] = ig;
      sparse["jg"] = jg;
      sparse["di"] = di;
      sparse["ggl"] = ggl;
      sparse["ggu"] = ggu;

      sparse_out << sparse;
      sparse_out.close();
   }
   break;

   default:
      throw "What type is this?\n";
   }
}


// ============================================================================
std::unique_ptr<dvector> profile_matrix::dot(const dvector &vector)
{
   auto result = std::make_unique <std::vector <double>>(vector.size());

   for (uint32_t i = 0; i < result->size(); i++)
   {
      (*result)[i] = di[i] * vector[i];

      uint32_t l = ig[i + 1] - ig[i];
      uint32_t k = i - 1;

      for (uint32_t j = 0; j < l; j++)
      {
         uint32_t index = ig[i] + j - 1;

         (*result)[i] += ggl[index] * vector[k];
         (*result)[k] += ggu[index] * vector[i];
      }
   }

   return result;
}

void profile_matrix::to_dense()
{
   dense_matrix.resize(dim);

   for (auto &row : dense_matrix)
      row.resize(dim);

   for (uint32_t i = 0; i < dim; i++)
   {
      dense_matrix[i][i] = di[i];

      uint32_t prof = ig[i + 1] - ig[i];
      uint32_t j0 = i - prof;
      uint32_t pos = 0;

      for (uint32_t j = j0; j < j0 + prof; j++)
      {
         uint32_t ij = ig[i] + pos;
         dense_matrix[i][j] = ggl[ij];
         dense_matrix[j][i] = ggu[ij];
         pos++;
      }
   }
}

void profile_matrix::save(matrix_type type, std::string path)
{
   switch (type)
   {
   case matrix_type::DENSE: {
      std::ofstream dense_out(path + "dense.json");
      nlohmann::json dense{};

      dense["dense"] = dense_matrix;

      dense_out << dense;
      dense_out.close();

      dense_matrix.clear();
   }
   break;

   case matrix_type::PROFILE: {
      std::ofstream profile_out(path + "profile.json");
      nlohmann::json profile{};

      profile["ig"] = ig;
      profile["di"] = di;
      profile["ggl"] = ggl;
      profile["ggu"] = ggu;

      profile_out << profile;
      profile_out.close();
   }
   break;
   }
}