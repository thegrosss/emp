#pragma once
#include <fstream>
#include "json.hpp"

using json = nlohmann::json;

struct JsonConfig
{
   std::string mesh_type;

   std::vector<double> x_lines;
   std::vector<double> y_lines;

   std::vector<uint32_t> x_interval_splits;
   std::vector<uint32_t> y_interval_splits;

   std::vector<double> x_interval_k;
   std::vector<double> y_interval_k;

   std::uint32_t nesting;

   struct Area
   {
      std::vector<double> x_range;
      std::vector<double> y_range;
      double lambda;
      double gamma;
      std::vector<uint32_t> boundary_conditions;
   };
   std::vector<Area> areas;

   static JsonConfig read_from_file(const std::string &filename)
   {
      std::ifstream file(filename);
      if (!file.is_open())
         throw std::runtime_error("Cannot open file: " + filename);

      json j;
      file >> j;

      JsonConfig config;

      auto &mesh = j["mesh"];
      config.mesh_type = mesh["type"];
      config.x_lines = mesh["x_lines"].get<std::vector<double>>();
      config.y_lines = mesh["y_lines"].get<std::vector<double>>();
      config.x_interval_splits = mesh["x_interval_splits"].get<std::vector<uint32_t>>();
      config.y_interval_splits = mesh["y_interval_splits"].get<std::vector<uint32_t>>();
      config.x_interval_k = mesh["x_interval_k"].get<std::vector<double>>();
      config.y_interval_k = mesh["y_interval_k"].get<std::vector<double>>();
      config.nesting = mesh["nesting"];

      for (const auto area_json : j["areas"])
      {
         Area area;
         area.x_range = area_json["x_range"].get<std::vector<double>>();
         area.y_range = area_json["y_range"].get<std::vector<double>>();
         area.lambda = area_json["lambda"];
         area.gamma = area_json["gamma"];
         
         auto bc = area_json["boundary_conditions"];
         area.boundary_conditions = {
            bc["bottom"],
            bc["left"],
            bc["top"],
            bc["right"]
         };

         config.areas.push_back(area);
      }

      return config;
   }
};