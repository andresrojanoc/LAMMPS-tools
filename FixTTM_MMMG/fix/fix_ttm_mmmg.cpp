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

/* ----------------------------------------------------------------------
   Contributing authors: Paul Crozier (SNL)
                         Carolyn Phillips (University of Michigan)
                         Andres Rojano (Modernized MMMG)
------------------------------------------------------------------------- */

#include "fix_ttm_mmmg.h"
#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "potential_file_reader.h"
#include "random_mars.h"
#include "update.h"
#include "utils.h"
#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;

// OFFSET avoids outside-of-box atoms being rounded to grid pts incorrectly
// SHIFT = 0.0 assigns atoms to lower-left grid pt
// SHIFT = 0.5 assigns atoms to nearest grid pt
// use SHIFT = 0.0 for now since it allows fix ave/chunk
//   to spatially average consistent with the TTM grid

static constexpr int OFFSET = 16384;
static constexpr double SHIFT = 0.0;

/* ---------------------------------------------------------------------- */

FixTTMMMMG::FixTTMMMMG(LAMMPS *lmp, int narg, char **arg)
    : FixTTM(lmp, 13, arg) { // 13 is to pass only 13 arguments to the parent
                             // class (FixTTM), and avoid keyword arguments.
  // Intializing custom flags
  cutoff_active = false;
  offset_active = false;
  cetable_active = false;
  ketable_active = false;
  set_active = infile_active = false;
  tinit = 0.0;
  average_electronic_temperature = 0.0;

  // Parsing all keywords,  standard (set, infile, outfile) and custom (cutoff,
  // offset, cetab, ketab).
  int iarg = 13;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "set") == 0) {
      if (iarg + 2 > narg)
        error->all(FLERR, "Illegal fix ttm_mmmg command");
      tinit = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      set_active = true;
      if (tinit <= 0.0)
        error->all(FLERR, "Fix ttm_mmmg initial temperature must be > 0.0");
      iarg += 2;
    } else if (strcmp(arg[iarg], "infile") == 0) {
      if (iarg + 2 > narg)
        error->all(FLERR, "Illegal fix ttm_mmmg command");
      infile = arg[iarg + 1];
      infile_active = true;
      iarg += 2;
    } else if (strcmp(arg[iarg], "outfile") == 0) {
      if (iarg + 3 > narg)
        error->all(FLERR, "Illegal fix ttm_mmmg command");
      outevery = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      outfile = arg[iarg + 2];
      iarg += 3;
    } else if (strcmp(arg[iarg], "cutoff") == 0) {
      if (iarg + 1 > narg)
        error->all(FLERR, "Illegal fix ttm_mmmg command");
      cutoff_active = true;
      iarg += 1;
    } else if (strcmp(arg[iarg], "offset") == 0) {
      if (iarg + 2 > narg)
        error->all(FLERR, "Illegal fix ttm_mmmg command");
      offset_active = true;
      time_offset = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "cetab") == 0) {
      if (iarg + 2 > narg)
        error->all(FLERR, "Illegal fix ttm_mmmg command");
      cetable_active = true;
      TableInterpReader(arg[iarg + 1], "ce");
      iarg += 2;
    } else if (strcmp(arg[iarg], "ketab") == 0) {
      if (iarg + 2 > narg)
        error->all(FLERR, "Illegal fix ttm_mmmg command");
      ketable_active = true;
      TableInterpReader(arg[iarg + 1], "ke");
      iarg += 2;
    } else {
      error->all(FLERR, "Illegal fix ttm_mmmg command");
    }
  }

  // If the cetable is activated, set the electronic_specific_heat value.
  if (cetable_active && set_active) {
    electronic_specific_heat =
        LinearInterpolate(tinit, "ce"); // To check stability_criterion
  }

  // If the ketable is activated, set the electronic_thermal_conductivity value.
  if (ketable_active && set_active) {
    electronic_thermal_conductivity =
        LinearInterpolate(tinit, "ke"); // To check stability_criterion
  }

  // error check
  if (seed <= 0)
    error->all(FLERR, "Invalid random number seed in fix ttm_mmmg command");
  if (electronic_specific_heat <= 0.0)
    error->all(FLERR, "Fix ttm_mmmg electronic_specific_heat must be > 0.0");

  if (cetable_active) {
    for (int i = 0; i < ce_values.size(); i++) {
      if (ce_values[i] <= 0.0)
        error->all(FLERR, "Fix ttm_mmmg all electronic_specific_heat entries "
                          "in the file must be > 0.0");
    }
  }

  if (ketable_active) {
    for (int i = 0; i < ke_values.size(); i++) {
      if (ke_values[i] <= 0.0)
        error->all(FLERR, "Fix ttm_mmmg all electronic_thermal_conductivity "
                          "entries in the file must be > 0.0");
    }
  }

  if (offset_active) {
    if (time_offset < 0)
      error->all(FLERR, "Fix ttm_mmmg offset must be >= 0");
  }

  if (offset_active && cutoff_active)
    error->all(
        FLERR,
        "Fix ttm_mmmg cannot have cutoff and offset active at the same time");

  if (electronic_density <= 0.0)
    error->all(FLERR, "Fix ttm_mmmg electronic_density must be > 0.0");
  if (electronic_thermal_conductivity < 0.0)
    error->all(FLERR,
               "Fix ttm_mmmg electronic_thermal_conductivity must be >= 0.0");
  if (gamma_p <= 0.0)
    error->all(FLERR, "Fix ttm_mmmg gamma_p must be > 0.0");
  if (gamma_s < 0.0)
    error->all(FLERR, "Fix ttm_mmmg gamma_s must be >= 0.0");
  if (v_0 < 0.0)
    error->all(FLERR, "Fix ttm_mmmg v_0 must be >= 0.0");
  if (nxgrid <= 0 || nygrid <= 0 || nzgrid <= 0)
    error->all(FLERR, "Fix ttm_mmmg grid sizes must be > 0");
}

/* ---------------------------------------------------------------------- */

FixTTMMMMG::~FixTTMMMMG() {
  // Parent class handles RNG and grid deallocation.
}

/* ---------------------------------------------------------------------- */

void FixTTMMMMG::post_force(int /*vflag*/) {
  int ix, iy, iz;
  double gamma1, gamma2;

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double *boxlo = domain->boxlo;
  double dxinv = nxgrid / domain->xprd;
  double dyinv = nygrid / domain->yprd;
  double dzinv = nzgrid / domain->zprd;

  update->update_time();
  double current_time = update->atime; // This is to look for the current time.
  gamma_offset = (offset_active && (current_time < time_offset)) ? 0 : 1;

  /*
  // Debug: Evaluates current simulation time

  if (comm->me == 0) {
    utils::logmesg(
        lmp, fmt::format("DEBUG: Time = {}, time_offset {}, gamma_offset = {}
  \n", update->atime, time_offset, gamma_offset));
  }
  */

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      ix = static_cast<int>((x[i][0] - boxlo[0]) * dxinv + shift) - OFFSET;
      iy = static_cast<int>((x[i][1] - boxlo[1]) * dyinv + shift) - OFFSET;
      iz = static_cast<int>((x[i][2] - boxlo[2]) * dzinv + shift) - OFFSET;
      if (ix < 0)
        ix += nxgrid;
      if (iy < 0)
        iy += nygrid;
      if (iz < 0)
        iz += nzgrid;
      if (ix >= nxgrid)
        ix -= nxgrid;
      if (iy >= nygrid)
        iy -= nygrid;
      if (iz >= nzgrid)
        iz -= nzgrid;
      if (T_electron[iz][iy][ix] < 0)
        error->one(FLERR, "Electronic temperature dropped below zero");

      double tsqrt = sqrt(T_electron[iz][iy][ix]);
      gamma1 = gfactor1[type[i]];
      double vsq = v[i][0] * v[i][0] + v[i][1] * v[i][1] + v[i][2] * v[i][2];

      gamma_cutoff = 1.0;

      // Evaluates whether the cutoff approach applies in the computation.
      if (vsq > v_0_sq) {
        if (cutoff_active) {
          gamma1 *=
              (gamma_s / gamma_p); // Avoids gamma_p for fast-moving atoms.
          gamma_cutoff = 0;
        } else {
          gamma1 *= ((gamma_p * gamma_offset) + gamma_s) /
                    gamma_p; // Standard ttm approach.
        }
      } else {                  // vsq <= v_0_sq
        gamma1 *= gamma_offset; // This decouples gamma_p for slow-moving atoms
                                // when time_offset applies.
      }

      gamma2 = gfactor2[type[i]] * tsqrt * gamma_cutoff * gamma_offset;
      // Modifies gamma2 to disable stochastic forces for fast-moving atoms and
      // at simulation times before time_offset.

      flangevin[i][0] = gamma1 * v[i][0] + gamma2 * (random->uniform() - 0.5);
      flangevin[i][1] = gamma1 * v[i][1] + gamma2 * (random->uniform() - 0.5);
      flangevin[i][2] = gamma1 * v[i][2] + gamma2 * (random->uniform() - 0.5);

      f[i][0] += flangevin[i][0];
      f[i][1] += flangevin[i][1];
      f[i][2] += flangevin[i][2];
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixTTMMMMG::end_of_step() {
  int ix, iy, iz;
  double **x = atom->x;
  double **v = atom->v;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double *boxlo = domain->boxlo;
  double dxinv = nxgrid / domain->xprd;
  double dyinv = nygrid / domain->yprd;
  double dzinv = nzgrid / domain->zprd;

  for (iz = 0; iz < nzgrid; iz++)
    for (iy = 0; iy < nygrid; iy++)
      for (ix = 0; ix < nxgrid; ix++)
        net_energy_transfer[iz][iy][ix] = 0.0;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      ix = static_cast<int>((x[i][0] - boxlo[0]) * dxinv + shift) - OFFSET;
      iy = static_cast<int>((x[i][1] - boxlo[1]) * dyinv + shift) - OFFSET;
      iz = static_cast<int>((x[i][2] - boxlo[2]) * dzinv + shift) - OFFSET;
      if (ix < 0)
        ix += nxgrid;
      if (iy < 0)
        iy += nygrid;
      if (iz < 0)
        iz += nzgrid;
      if (ix >= nxgrid)
        ix -= nxgrid;
      if (iy >= nygrid)
        iy -= nygrid;
      if (iz >= nzgrid)
        iz -= nzgrid;

      net_energy_transfer[iz][iy][ix] +=
          (flangevin[i][0] * v[i][0] + flangevin[i][1] * v[i][1] +
           flangevin[i][2] * v[i][2]);
    }

  outflag = 0;
  MPI_Allreduce(&net_energy_transfer[0][0][0],
                &net_energy_transfer_all[0][0][0], ngridtotal, MPI_DOUBLE,
                MPI_SUM, world);

  double dx = domain->xprd / nxgrid;
  double dy = domain->yprd / nygrid;
  double dz = domain->zprd / nzgrid;
  double del_vol = dx * dy * dz;

  // num_inner_timesteps = # of inner steps (thermal solves)
  // required this MD step to maintain a stable explicit solve

  int num_inner_timesteps = 1;
  double inner_dt = update->dt;

  double stability_criterion =
      1.0 - 2.0 * inner_dt / (electronic_specific_heat * electronic_density) *
                (electronic_thermal_conductivity *
                 (1.0 / dx / dx + 1.0 / dy / dy + 1.0 / dz / dz));

  if (stability_criterion < 0.0) {
    inner_dt = 0.5 * (electronic_specific_heat * electronic_density) /
               (electronic_thermal_conductivity *
                (1.0 / dx / dx + 1.0 / dy / dy + 1.0 / dz / dz));
    num_inner_timesteps = static_cast<int>(update->dt / inner_dt) + 1;
    inner_dt = update->dt / double(num_inner_timesteps);
    if (num_inner_timesteps > 1000000)
      error->warning(FLERR, "Too many inner timesteps in fix ttm_mmmg");
  }

  // finite difference iterations to update T_electron

  for (int istep = 0; istep < num_inner_timesteps; istep++) {

    for (iz = 0; iz < nzgrid; iz++)
      for (iy = 0; iy < nygrid; iy++)
        for (ix = 0; ix < nxgrid; ix++)
          T_electron_old[iz][iy][ix] = T_electron[iz][iy][ix];

    // compute new electron T profile

    for (iz = 0; iz < nzgrid; iz++)
      for (iy = 0; iy < nygrid; iy++)
        for (ix = 0; ix < nxgrid; ix++) {
          int xright = ix + 1;
          int yright = iy + 1;
          int zright = iz + 1;
          if (xright == nxgrid)
            xright = 0;
          if (yright == nygrid)
            yright = 0;
          if (zright == nzgrid)
            zright = 0;
          int xleft = ix - 1;
          int yleft = iy - 1;
          int zleft = iz - 1;
          if (xleft == -1)
            xleft = nxgrid - 1;
          if (yleft == -1)
            yleft = nygrid - 1;
          if (zleft == -1)
            zleft = nzgrid - 1;

          variable_electronic_specific_heat =
              cetable_active ? LinearInterpolate(T_electron[iz][iy][ix], "ce")
                             : electronic_specific_heat;

          variable_electronic_thermal_conductivity =
              ketable_active ? LinearInterpolate(T_electron[iz][iy][ix], "ke")
                             : electronic_thermal_conductivity;

          /*
          // Debug: Evaluates the interpolation of the
          variable_electronic_specific_heat:

            if (comm->me == 0 && update->ntimestep % 10 == 0 && ix == 0 &&
                iy == 0 && iz == 0) {
              utils::logmesg(
                  lmp, fmt::format("DEBUG: Step {} | T={} K | ce={} | ke={}\n",
                                   update->ntimestep, T_electron[0][0][0],
                                   variable_electronic_specific_heat,
                                   variable_electronic_thermal_conductivity));
            }
          */

          T_electron[iz][iy][ix] =
              T_electron_old[iz][iy][ix] +
              inner_dt /
                  (variable_electronic_specific_heat * electronic_density) *
                  (variable_electronic_thermal_conductivity *

                       ((T_electron_old[iz][iy][xright] +
                         T_electron_old[iz][iy][xleft] -
                         2.0 * T_electron_old[iz][iy][ix]) /
                            dx / dx +
                        (T_electron_old[iz][yright][ix] +
                         T_electron_old[iz][yleft][ix] -
                         2.0 * T_electron_old[iz][iy][ix]) /
                            dy / dy +
                        (T_electron_old[zright][iy][ix] +
                         T_electron_old[zleft][iy][ix] -
                         2.0 * T_electron_old[iz][iy][ix]) /
                            dz / dz) -

                   (net_energy_transfer_all[iz][iy][ix]) / del_vol);
        }
  }

  // output of grid electron temperatures to file

  if (!outfile.empty() && (update->ntimestep % outevery == 0))
    write_electron_temperatures(
        fmt::format("{}.{}", outfile, update->ntimestep));
}

/* ---------------------------------------------------------------------- */

void FixTTMMMMG::read_electron_temperatures(const std::string &filename) {
  if (comm->me == 0) {

    int ***T_initial_set;
    memory->create(T_initial_set, nzgrid, nygrid, nxgrid,
                   "ttm_mmmg:T_initial_set");
    memset(&T_initial_set[0][0][0], 0, ngridtotal * sizeof(int));

    // read initial electron temperature values from file
    bigint nread = 0;
    double sum = 0.0;

    try {
      PotentialFileReader reader(lmp, filename, "electron temperature grid");

      while (nread < ngridtotal) {
        // reader will skip over comment-only lines
        auto values = reader.next_values(4);
        ++nread;

        int ix = values.next_int() - 1;
        int iy = values.next_int() - 1;
        int iz = values.next_int() - 1;
        double T_tmp = values.next_double();

        // check correctness of input data

        if ((ix < 0) || (ix >= nxgrid) || (iy < 0) || (iy >= nygrid) ||
            (iz < 0) || (iz >= nzgrid))
          throw TokenizerException(
              "Fix ttm_mmmg invalid grid index in fix ttm grid file", "");

        if (T_tmp < 0.0)
          throw TokenizerException(
              "Fix ttm_mmmg electron temperatures must be > 0.0", "");

        T_electron[iz][iy][ix] = T_tmp;
        T_initial_set[iz][iy][ix] = 1;
        sum += T_tmp;
      }
    } catch (std::exception &e) {
      error->one(FLERR, e.what());
    }

    // check completeness of input data

    for (int iz = 0; iz < nzgrid; iz++)
      for (int iy = 0; iy < nygrid; iy++)
        for (int ix = 0; ix < nxgrid; ix++)
          if (T_initial_set[iz][iy][ix] == 0)
            error->all(FLERR,
                       "Fix ttm_mmmg infile did not set all temperatures");

    memory->destroy(T_initial_set);
    average_electronic_temperature = sum / ngridtotal;

    if (cetable_active)
      electronic_specific_heat = LinearInterpolate(
          average_electronic_temperature, "ce"); // To check stability_criterion

    if (ketable_active)
      electronic_thermal_conductivity = LinearInterpolate(
          average_electronic_temperature, "ke"); // To check stability_criterion
  }

  MPI_Bcast(&T_electron[0][0][0], ngridtotal, MPI_DOUBLE, 0, world);
  MPI_Bcast(&average_electronic_temperature, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&electronic_specific_heat, 1, MPI_DOUBLE, 0, world);
}

/* ---------------------------------------------------------------------- */

double FixTTMMMMG::compute_vector(int n) {
  if (outflag == 0) {
    e_energy = 0.0;
    transfer_energy = 0.0;

    double dx = domain->xprd / nxgrid;
    double dy = domain->yprd / nygrid;
    double dz = domain->zprd / nzgrid;
    double del_vol = dx * dy * dz;

    for (int iz = 0; iz < nzgrid; iz++) {
      for (int iy = 0; iy < nygrid; iy++) {
        for (int ix = 0; ix < nxgrid; ix++) {
          variable_electronic_specific_heat =
              cetable_active ? LinearInterpolate(T_electron[iz][iy][ix], "ce")
                             : electronic_specific_heat;
          e_energy += T_electron[iz][iy][ix] *
                      variable_electronic_specific_heat * electronic_density *
                      del_vol;
          transfer_energy += net_energy_transfer_all[iz][iy][ix] * update->dt;
        }
      }
    }
    outflag = 1;
  }

  if (n == 0)
    return e_energy;
  if (n == 1)
    return transfer_energy;
  return 0.0;
}

/* ---------------------------------------------------------------------- */

void FixTTMMMMG::TableInterpReader(const std::string &filename,
                                   const std::string &keyword) {

  std::vector<double> &temp_vals =
      (keyword == "ce") ? temp_ce_values : temp_ke_values;
  std::vector<double> &dtemp_vals =
      (keyword == "ce") ? dtemp_ce_values : dtemp_ke_values;
  std::vector<double> &y_vals = (keyword == "ce") ? ce_values : ke_values;
  std::vector<double> &dy_vals = (keyword == "ce") ? dce_values : dke_values;

  if (comm->me == 0) {
    std::string table_label = (keyword == "ce") ? "specific heat table"
                                                : "thermal conductivity table";
    PotentialFileReader reader(lmp, filename, table_label);
    while (char *line = reader.next_line()) {
      double temp_value, y_value;
      if (sscanf(line, "%lg %lg", &temp_value, &y_value) == 2) {
        temp_vals.push_back(temp_value);
        y_vals.push_back(y_value);
      }
    }

    // Pre-calculate deltas
    int nsize_table = static_cast<int>(temp_vals.size());
    dtemp_vals.resize(nsize_table);
    dy_vals.resize(nsize_table);
    for (int i = 0; i < nsize_table - 1; i++) {
      dtemp_vals[i] = temp_vals[i + 1] - temp_vals[i];
      dy_vals[i] = y_vals[i + 1] - y_vals[i];
    }
  }

  int nsize_table = static_cast<int>(temp_vals.size());
  MPI_Bcast(&nsize_table, 1, MPI_INT, 0, world);

  if (comm->me != 0) {
    temp_vals.resize(nsize_table);
    y_vals.resize(nsize_table);
    dtemp_vals.resize(nsize_table);
    dy_vals.resize(nsize_table);
  }

  MPI_Bcast(temp_vals.data(), nsize_table, MPI_DOUBLE, 0, world);
  MPI_Bcast(y_vals.data(), nsize_table, MPI_DOUBLE, 0, world);
  MPI_Bcast(dtemp_vals.data(), nsize_table, MPI_DOUBLE, 0, world);
  MPI_Bcast(dy_vals.data(), nsize_table, MPI_DOUBLE, 0, world);
}

/* ---------------------------------------------------------------------- */

double FixTTMMMMG::LinearInterpolate(double temp, const std::string &keyword) {

  const std::vector<double> &temp_vals =
      (keyword == "ce") ? temp_ce_values : temp_ke_values;
  const std::vector<double> &dtemp_vals =
      (keyword == "ce") ? dtemp_ce_values : dtemp_ke_values;
  const std::vector<double> &y_vals = (keyword == "ce") ? ce_values : ke_values;
  const std::vector<double> &dy_vals =
      (keyword == "ce") ? dce_values : dke_values;

  if (temp_vals.empty())
    error->one(FLERR, "The table {} is empty", keyword);

  int lo = 0;
  int hi = static_cast<int>(temp_vals.size()) - 1;

  if (temp <= temp_vals[lo])
    // Exact match or lower than lower bound
    return y_vals[lo];
  if (temp >= temp_vals[hi])
    // Exact match or higher than upper bound
    return y_vals[hi];

  while (hi - lo > 1) {
    int mid = (lo + hi) / 2;
    if (temp < temp_vals[mid])
      hi = mid;
    else
      lo = mid;
  }

  // Use pre-calculated deltas
  return y_vals[lo] + (temp - temp_vals[lo]) * (dy_vals[lo] / dtemp_vals[lo]);
}
