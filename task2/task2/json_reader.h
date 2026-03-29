#pragma once
#include <fstream>
#include "json.hpp"

using json = nlohmann::json;

struct JsonConfig
{
   std::string mfe_method;
   bool use_relax;

   struct Area
   {
      std::string type;

      double start;
      double end;
      double sigma;

      uint32_t nx;
      double kx;

      int nesting;
   };
   Area area;

   struct Time
   {
      std::string type;
      double start;
      double end;
      uint32_t nt;
      double kt;
      int nesting;
   };
   Time time;

   static JsonConfig from_file(const std::string& path)
   {
      std::ifstream in(path);
      if (!in.is_open())
         throw "Can't open file";
      
      json j;
      in >> j;

      JsonConfig config;

      config.mfe_method = j["mfe_method"];
      config.use_relax = j["relaxation"];

      config.area.type = j["area"]["type"];
      config.area.start = j["area"]["limits"][0];
      config.area.end = j["area"]["limits"][1];
      config.area.sigma = j["area"]["sigma"];
      config.area.nx = j["area"]["elements_count"];
      config.area.kx = j["area"]["discharge_coefficient"];
      config.area.nesting = j["area"]["nesting"];

      config.time.type = j["time"]["type"];
      config.time.start = j["time"]["limits"][0];
      config.time.end = j["time"]["limits"][1];
      config.time.nt = j["time"]["elements_count"];
      config.time.kt = j["time"]["discharge_coefficient"];
      config.time.nesting = j["time"]["nesting"];

      return config;
   }
};