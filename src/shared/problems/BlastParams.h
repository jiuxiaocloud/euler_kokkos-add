#ifndef BLAST_PARAMS_H_
#define BLAST_PARAMS_H_

#include "utils/config/ConfigMap.h"

struct BlastParams
{

  // blast problem parameters
  real_t blast_radius;
  real_t blast_center_x;
  real_t blast_center_y;
  real_t blast_center_z;
  real_t blast_density_in;
  real_t blast_density_out;
  real_t blast_pressure_in;
  real_t blast_pressure_out;
#if COP
  // NOTE: add compoent problem parameters
  bool if_cop;
  real_t cop_radius;
  real_t cop_center_x;
  real_t cop_center_y;
  real_t cop_pressure_in;
  real_t cop_pressure_out;
  real_t cop_y1_in;
  real_t cop_y1_out;
  real_t cop_density_in;
  real_t cop_density_out;
#endif
  BlastParams(ConfigMap &configMap)
  {

    double xmin = configMap.getFloat("mesh", "xmin", 0.0);
    double ymin = configMap.getFloat("mesh", "ymin", 0.0);
    double zmin = configMap.getFloat("mesh", "zmin", 0.0);

    double xmax = configMap.getFloat("mesh", "xmax", 1.0);
    double ymax = configMap.getFloat("mesh", "ymax", 1.0);
    double zmax = configMap.getFloat("mesh", "zmax", 1.0);

    blast_radius = configMap.getFloat("blast", "blast_radius", (xmin + xmax) / 2.0 / 10);
    blast_center_x = configMap.getFloat("blast", "blast_center_x", (xmin + xmax) / 2);
    blast_center_y = configMap.getFloat("blast", "blast_center_y", (ymin + ymax) / 2);
    blast_center_z = configMap.getFloat("blast", "blast_center_z", (zmin + zmax) / 2);
    blast_density_in = configMap.getFloat("blast", "blast_density_in", 1.0);
    blast_density_out = configMap.getFloat("blast", "blast_density_out", 1.2);
    blast_pressure_in = configMap.getFloat("blast", "blast_pressure_in", 10.0);
    blast_pressure_out = configMap.getFloat("blast", "blast_pressure_out", 0.1);
#if COP
    // NOTE: add compoent read function
    if_cop = configMap.getFloat("blast", "if_cop", 0);
    cop_radius = configMap.getFloat("blast", "cop_radius", (xmin + xmax) / 2.0 / 10);
    cop_center_x = configMap.getFloat("blast", "cop_center_x", (xmin + xmax) / 2);
    cop_center_y = configMap.getFloat("blast", "cop_center_y", (ymin + ymax) / 2);
    cop_density_in = configMap.getFloat("blast", "cop_density_in", 1.0);
    cop_density_out = configMap.getFloat("blast", "cop_density_out", 1.2);
    cop_pressure_in = configMap.getFloat("blast", "cop_pressure_in", 10.0);
    cop_pressure_out = configMap.getFloat("blast", "cop_pressure_out", 0.1);
    cop_y1_in = configMap.getFloat("blast", "cop_y1_in", 0.0);
    cop_y1_out = configMap.getFloat("blast", "cop_y1_out", 1.0);
    cop_density_in = configMap.getFloat("blast", "cop_density_in", 1.2);
#endif
  }

}; // struct BlastParams

#endif // BLAST_PARAMS_H_
