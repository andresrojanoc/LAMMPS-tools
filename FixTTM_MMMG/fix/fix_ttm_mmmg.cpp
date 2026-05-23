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
  cetable_active = false;
  set_active = infile_active = false;
  tinit = 0.0;
  average_electronic_temperature = 0.0;

  // Parsing all keywords,  standard (set, infile, outfile) and custom (cutoff,
  // table).
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
    } else if (strcmp(arg[iarg], "table") == 0) {
      if (iarg + 2 > narg)
        error->all(FLERR, "Illegal fix ttm_mmmg command");
      cetable_active = true;
      TableInterpReader(arg[iarg + 1], "Ce");
      iarg += 2;
    } else {
      error->all(FLERR, "Illegal fix ttm_mmmg command");
    }
  }

  // If the cetable is activated, set the electronic_specific_heat value.
  if (cetable_active && set_active) {
    electronic_specific_heat =
        LinearInterpolate(tinit, "Ce"); // To check stability_criterion
  }

  // error check
  if (seed <= 0)
    error->all(FLERR, "Invalid random number seed in fix ttm_mmmg command");
  if (electronic_specific_heat <= 0.0)
    error->all(FLERR, "Fix ttm_mmmg electronic_specific_heat must be > 0.0");

  if (cetable_active) {
    for (int i = 0; i < Ce_values.size(); i++) {
      if (Ce_values[i] <= 0.0)
        error->all(FLERR, "Fix ttm_mmmg all electronic_specific_heat entries "
                          "in the file must be > 0.0");
    }
  }

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

      g_unit = 1.0;

      // Evaluates whether the cutoff approach applies in the computation.
      if (cutoff_active) {
        if (vsq > v_0_sq) {
          gamma1 *=
              (gamma_s / gamma_p); // Avoids gamma_p for fast-moving atoms.
          g_unit = 0;
        }
      } else {
        if (vsq > v_0_sq)
          gamma1 *= (gamma_p + gamma_s) / gamma_p; // Standard ttm approach.
      }

      gamma2 = gfactor2[type[i]] * tsqrt *
               g_unit; // Modifies gamma2 to disable stochastic forces for
                       // fast-moving atoms.

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
              cetable_active ? LinearInterpolate(T_electron[iz][iy][ix], "Ce")
                             : electronic_specific_heat;

          /*
          //Debug: Evaluates the interpolation of the
          variable_electronic_specific_heat:

          if (comm->me == 0 &&
          update->ntimestep % 10 == 0 && ix == 0 && iy == 0 && iz == 0) {
            utils::logmesg(lmp,
                           fmt::format("DEBUG: Step {} | T={} K | Ce={}\n",
                                       update->ntimestep, T_electron[0][0][0],
                                       variable_electronic_specific_heat));
          }
          */

          T_electron[iz][iy][ix] =
              T_electron_old[iz][iy][ix] +
              inner_dt /
                  (variable_electronic_specific_heat * electronic_density) *
                  (electronic_thermal_conductivity *

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
          average_electronic_temperature, "Ce"); // To check stability_criterion
  }

  MPI_Bcast(&T_electron[0][0][0], ngridtotal, MPI_DOUBLE, 0, world);
  MPI_Bcast(&average_electronic_temperature, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&electronic_specific_heat, 1, MPI_DOUBLE, 0, world);
}

/* ---------------------------------------------------------------------- */

void FixTTMMMMG::TableInterpReader(const std::string &filename,
                                   const std::string &keyword) {
  if (comm->me == 0) {
    PotentialFileReader reader(lmp, filename, "specific heat table");
    while (char *line = reader.next_line()) {
      double Temp_Ce_value, Ce_value;
      if (sscanf(line, "%lg %lg", &Temp_Ce_value, &Ce_value) == 2) {
        Temp_Ce_values.push_back(Temp_Ce_value);
        Ce_values.push_back(Ce_value);
      }
    }

    // Pre-calculate deltas (Recovered Article Logic)
    int nsize_cetable = static_cast<int>(Temp_Ce_values.size());
    dTemp_Ce_values.resize(nsize_cetable);
    dCe_values.resize(nsize_cetable);
    for (int i = 0; i < nsize_cetable - 1; i++) {
      dTemp_Ce_values[i] = Temp_Ce_values[i + 1] - Temp_Ce_values[i];
      dCe_values[i] = Ce_values[i + 1] - Ce_values[i];
    }
  }

  int nsize_cetable = static_cast<int>(Temp_Ce_values.size());
  MPI_Bcast(&nsize_cetable, 1, MPI_INT, 0, world);

  if (comm->me != 0) {
    Temp_Ce_values.resize(nsize_cetable);
    Ce_values.resize(nsize_cetable);
    dTemp_Ce_values.resize(nsize_cetable);
    dCe_values.resize(nsize_cetable);
  }
  MPI_Bcast(Temp_Ce_values.data(), nsize_cetable, MPI_DOUBLE, 0, world);
  MPI_Bcast(Ce_values.data(), nsize_cetable, MPI_DOUBLE, 0, world);
  MPI_Bcast(dTemp_Ce_values.data(), nsize_cetable, MPI_DOUBLE, 0, world);
  MPI_Bcast(dCe_values.data(), nsize_cetable, MPI_DOUBLE, 0, world);
}

/* ---------------------------------------------------------------------- */

double FixTTMMMMG::LinearInterpolate(double Temp_Ce,
                                     const std::string &keyword) {
  if (Temp_Ce_values.empty())
    error->one(FLERR, "The table {} is empty", keyword);

  int lo = 0;
  int hi = static_cast<int>(Temp_Ce_values.size()) - 1;

  if (Temp_Ce <= Temp_Ce_values[lo])
    // Exact match or lower than lower bound
    return Ce_values[lo];
  if (Temp_Ce >= Temp_Ce_values[hi])
    // Exact match or higher than upper bound
    return Ce_values[hi];

  while (hi - lo > 1) {
    int mid = (lo + hi) / 2;
    if (Temp_Ce < Temp_Ce_values[mid])
      hi = mid;
    else
      lo = mid;
  }

  // Use pre-calculated deltas (Recovered Article Logic)
  return Ce_values[lo] + (Temp_Ce - Temp_Ce_values[lo]) *
                             (dCe_values[lo] / dTemp_Ce_values[lo]);
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
              cetable_active ? LinearInterpolate(T_electron[iz][iy][ix], "Ce")
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