/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS
// clang-format off
FixStyle(ttm_mmmg,FixTTMMMMG);
// clang-format on
#else

#ifndef LMP_FIX_TTM_MMMG_H
#define LMP_FIX_TTM_MMMG_H

#include "fix_ttm.h"
#include <string>
#include <vector>

namespace LAMMPS_NS {

class FixTTMMMMG : public FixTTM {
public:
  FixTTMMMMG(class LAMMPS *, int, char **);
  ~FixTTMMMMG() override;

  void post_force(int) override;
  void end_of_step() override;
  double compute_vector(int) override;

protected:
  // --- MMMG Specific Flags (bool) ---
  bool cutoff_active, offset_active, cetable_active, ketable_active;
  bool set_active, infile_active;

  double time_offset;
  double gamma_cutoff, gamma_offset;
  double average_electronic_temperature;
  double variable_electronic_specific_heat,
      variable_electronic_thermal_conductivity;

  void read_electron_temperatures(const std::string &filename) override;

  // Tabular Specific Heat data
  std::vector<double> temp_ce_values;
  std::vector<double> ce_values;
  std::vector<double> dtemp_ce_values;
  std::vector<double> dce_values;

  // Tabular Thermal Conductivity data
  std::vector<double> temp_ke_values;
  std::vector<double> ke_values;
  std::vector<double> dtemp_ke_values;
  std::vector<double> dke_values;

  // Internal Logic Methods
  void TableInterpReader(const std::string &filename,
                         const std::string &keyword);
  double LinearInterpolate(double temp, const std::string &keyword);
};

} // namespace LAMMPS_NS

#endif
#endif
