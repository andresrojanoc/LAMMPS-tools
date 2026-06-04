// clang-format off
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
------------------------------------------------------------------------- */

#include "fix_ttm_cascade.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "grid3d.h"
#include "memory.h"
#include "neighbor.h"
#include "potential_file_reader.h"
#include "random_mars.h"
#include "safe_pointers.h"
#include "tokenizer.h"
#include "update.h"

#include <cmath>
#include <cstring>
#include <exception>

using namespace LAMMPS_NS;
using namespace FixConst;

static constexpr int MAXLINE = 256;
static constexpr int CHUNK = 1024;

// OFFSET avoids outside-of-box atoms being rounded to grid pts incorrectly

static constexpr int OFFSET = 16384;

/* ---------------------------------------------------------------------- */

FixTTMCascade::FixTTMCascade(LAMMPS *lmp, int narg, char **arg)
    : FixTTM(lmp, 13,
             arg) // 13 is to pass only 13 arguments to the parent class
                  // (FixTTMGrid and FixTTM), and avoid keyword arguments.
{
  cutoff_active = false;
  offset_active = false;
  cetable_active = false;
  ketable_active = false;
  set_active = infile_active = false;
  tinit = 0.0;
  average_electronic_temperature = 0.0;

  pergrid_flag = 1;
  pergrid_freq = 1;
  restart_file = 1;

  // Parsing all keywords,  standard (set and infile) and custom (cutoff,
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
      tableinterpreader(arg[iarg + 1], "ce");
      iarg += 2;
    } else if (strcmp(arg[iarg], "ketab") == 0) {
      if (iarg + 2 > narg)
        error->all(FLERR, "Illegal fix ttm_mmmg command");
      ketable_active = true;
      tableinterpreader(arg[iarg + 1], "ke");
      iarg += 2;
    } else {
      error->all(FLERR, "Illegal fix ttm_mmmg command");
    }
  }

  // If the cetable is activated, set the electronic_specific_heat value.
  if (cetable_active && set_active) {
    electronic_specific_heat =
        linearinterpolation(tinit, "ce"); // To check stability_criterion
  }

  // If the ketable is activated, set the electronic_thermal_conductivity value.
  if (ketable_active && set_active) {
    electronic_thermal_conductivity =
        linearinterpolation(tinit, "ke"); // To check stability_criterion
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

  if (outfile.size() > 0)
    error->all(FLERR, Error::NOPOINTER,
               "Fix ttm/cascade does not support outfile option - "
               "use dump grid command or restart files instead");

  skin_original = neighbor->skin;
}

/* ---------------------------------------------------------------------- */

FixTTMCascade::~FixTTMCascade()
{
  FixTTMCascade::deallocate_grid();
  deallocate_flag = 1;
}

/* ---------------------------------------------------------------------- */

void FixTTMCascade::post_constructor()
{
  // allocate global grid on each proc
  // needs to be done in post_contructor() beccause is virtual method

  allocate_grid();

  // initialize electron temperatures on grid

  int ix,iy,iz;
  for (iz = nzlo_out; iz <= nzhi_out; iz++)
    for (iy = nylo_out; iy <= nyhi_out; iy++)
      for (ix = nxlo_out; ix <= nxhi_out; ix++)
        T_electron[iz][iy][ix] = tinit;

  // zero net_energy_transfer
  // in case compute_vector accesses it on timestep 0

  outflag = 0;
  memset(&net_energy_transfer[nzlo_out][nylo_out][nxlo_out],0,
         ngridout*sizeof(double));

  // set initial electron temperatures from user input file
  // communicate new T_electron values to ghost grid points

  if (!infile.empty()) {
    read_electron_temperatures(infile);
    grid->forward_comm(Grid3d::FIX,this,0,1,sizeof(double), grid_buf1,grid_buf2,MPI_DOUBLE);
  }
}

/* ---------------------------------------------------------------------- */

void FixTTMCascade::init()
{
  FixTTM::init();

  if (neighbor->skin > skin_original)
    error->all(FLERR, Error::NOLASTLINE,
               "Cannot extend neighbor skin after fix ttm/cascade defined");
}

/* ---------------------------------------------------------------------- */

void FixTTMCascade::post_force(int /*vflag*/)
{
  int ix,iy,iz;
  double gamma1,gamma2;

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double *boxlo = domain->boxlo;
  double dxinv = nxgrid/domain->xprd;
  double dyinv = nygrid/domain->yprd;
  double dzinv = nzgrid/domain->zprd;

  // apply damping and thermostat to all atoms in fix group

  int flag = 0;

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
      ix = static_cast<int> ((x[i][0]-boxlo[0])*dxinv + OFFSET) - OFFSET;
      iy = static_cast<int> ((x[i][1]-boxlo[1])*dyinv + OFFSET) - OFFSET;
      iz = static_cast<int> ((x[i][2]-boxlo[2])*dzinv + OFFSET) - OFFSET;

      // flag if ix,iy,iz is not within my ghost cell range

      if (ix < nxlo_out || ix > nxhi_out || iy < nylo_out || iy > nyhi_out ||
          iz < nzlo_out || iz > nzhi_out) {
        flag = 1;
        continue;
      }

      if (T_electron[iz][iy][ix] < 0)
        error->one(FLERR,"Electronic temperature dropped below zero");

      double tsqrt = sqrt(T_electron[iz][iy][ix]);
      gamma1 = gfactor1[type[i]];
      double vsq = v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2];

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

      flangevin[i][0] = gamma1*v[i][0] + gamma2*(random->uniform()-0.5);
      flangevin[i][1] = gamma1*v[i][1] + gamma2*(random->uniform()-0.5);
      flangevin[i][2] = gamma1*v[i][2] + gamma2*(random->uniform()-0.5);

      f[i][0] += flangevin[i][0];
      f[i][1] += flangevin[i][1];
      f[i][2] += flangevin[i][2];
    }
  }

  if (flag) error->one(FLERR,"Out of range fix ttm/cascade atoms");
}

/* ---------------------------------------------------------------------- */

void FixTTMCascade::end_of_step(){
  int ix,iy,iz;

  double **x = atom->x;
  double **v = atom->v;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double *boxlo = domain->boxlo;
  double dxinv = nxgrid/domain->xprd;
  double dyinv = nygrid/domain->yprd;
  double dzinv = nzgrid/domain->zprd;
  double volgrid = 1.0 / (dxinv*dyinv*dzinv);

  double variable_electronic_specific_heat;
  double variable_electronic_thermal_conductivity;
  double el_th_diffusivity;
  double el_th_diffusivity_max = 0.0;
  double el_th_diffusivity_global_max = 0.0;

  outflag = 0;
  memset(&net_energy_transfer[nzlo_out][nylo_out][nxlo_out],0,
         ngridout*sizeof(double));

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      ix = static_cast<int> ((x[i][0]-boxlo[0])*dxinv + OFFSET) - OFFSET;
      iy = static_cast<int> ((x[i][1]-boxlo[1])*dyinv + OFFSET) - OFFSET;
      iz = static_cast<int> ((x[i][2]-boxlo[2])*dzinv + OFFSET) - OFFSET;

      net_energy_transfer[iz][iy][ix] +=
          (flangevin[i][0]*v[i][0] + flangevin[i][1]*v[i][1] +
           flangevin[i][2]*v[i][2]);
    }

  grid->reverse_comm(Grid3d::FIX,this,0,1,sizeof(double),
                     grid_buf1,grid_buf2,MPI_DOUBLE);

  // clang-format off

  // num_inner_timesteps = # of inner steps (thermal solves)
  // required this MD step to maintain a stable explicit solve

  int num_inner_timesteps = 1;
  double inner_dt = update->dt;

  if (cetable_active || ketable_active) {

    for (iz = nzlo_in; iz <= nzhi_in; iz++)
      for (iy = nylo_in; iy <= nyhi_in; iy++)
        for (ix = nxlo_in; ix <= nxhi_in; ix++) {
          variable_electronic_specific_heat = cetable_active ? linearinterpolation(T_electron[iz][iy][ix], "ce") : electronic_specific_heat;
          variable_electronic_thermal_conductivity = ketable_active ? linearinterpolation(T_electron[iz][iy][ix], "ke") : electronic_thermal_conductivity;
          el_th_diffusivity = variable_electronic_thermal_conductivity/(variable_electronic_specific_heat*electronic_density);
          if (el_th_diffusivity > el_th_diffusivity_max) el_th_diffusivity_max = el_th_diffusivity;
        }
    MPI_Allreduce(&el_th_diffusivity_max, &el_th_diffusivity_global_max, 1, MPI_DOUBLE, MPI_MAX, world);    
  }
  

  else {
    el_th_diffusivity_global_max = electronic_thermal_conductivity/(electronic_specific_heat*electronic_density);
  }

  double stability_criterion = 1.0 -
    2.0*inner_dt* el_th_diffusivity_global_max * (dxinv*dxinv + dyinv*dyinv + dzinv*dzinv);

  if (stability_criterion < 0.0) {
    inner_dt = 0.5/ el_th_diffusivity_global_max / (dxinv*dxinv + dyinv*dyinv + dzinv*dzinv);
    num_inner_timesteps = static_cast<int>(update->dt/inner_dt) + 1;
    inner_dt = update->dt/double(num_inner_timesteps);
    if (num_inner_timesteps > 1000000)
      error->warning(FLERR,"Too many inner timesteps in fix ttm/cascade");
  }


  // finite difference iterations to update T_electron

  for (int istep = 0; istep < num_inner_timesteps; istep++) {

    memcpy(&T_electron_old[nzlo_out][nylo_out][nxlo_out],
           &T_electron[nzlo_out][nylo_out][nxlo_out],ngridout*sizeof(double));

  // store thermal conductivity in a grid for rapid access           

  if(ketable_active){

    memset(&conductivity_xface[nzlo_out][nylo_out][nxlo_out],0, ngridout*sizeof(double));
    memset(&conductivity_yface[nzlo_out][nylo_out][nxlo_out],0, ngridout*sizeof(double));
    memset(&conductivity_zface[nzlo_out][nylo_out][nxlo_out],0, ngridout*sizeof(double));

    // x faces

    for (iz = nzlo_out; iz <= nzhi_out; iz++)
      for (iy = nylo_out; iy <= nyhi_out; iy++)
        for (ix = nxlo_out; ix < nxhi_out; ix++){
          conductivity_xface[iz][iy][ix] = linearinterpolation((T_electron_old[iz][iy][ix]+T_electron_old[iz][iy][ix+1])/2.0, "ke");
          if(conductivity_xface[iz][iy][ix]<=0) 
            error->all(FLERR,"Fix ttm/cascade: invalid conductivity at x-face ({},{},{}) value={}",ix, iy, iz,conductivity_xface[iz][iy][ix]);
        }

    // y faces

    for (iz = nzlo_out; iz <= nzhi_out; iz++)
      for (iy = nylo_out; iy < nyhi_out; iy++)
        for (ix = nxlo_out; ix <= nxhi_out; ix++){
          conductivity_yface[iz][iy][ix] = linearinterpolation((T_electron_old[iz][iy][ix]+T_electron_old[iz][iy+1][ix])/2.0, "ke");
          if(conductivity_yface[iz][iy][ix]<=0) 
            error->all(FLERR,"Fix ttm/cascade: invalid conductivity at y-face ({},{},{}) value={}",ix, iy, iz,conductivity_yface[iz][iy][ix]);
        }

    // z faces

    for (iz = nzlo_out; iz < nzhi_out; iz++)
      for (iy = nylo_out; iy <= nyhi_out; iy++)
        for (ix = nxlo_out; ix <= nxhi_out; ix++){
          conductivity_zface[iz][iy][ix] = linearinterpolation((T_electron_old[iz][iy][ix]+T_electron_old[iz+1][iy][ix])/2.0, "ke");
          if(conductivity_zface[iz][iy][ix]<=0) 
            error->all(FLERR,"Fix ttm/cascade: invalid conductivity at z-face ({},{},{}) value={}",ix, iy, iz,conductivity_zface[iz][iy][ix]);
        }
  }

    // compute new electron T profile

    for (iz = nzlo_in; iz <= nzhi_in; iz++)
      for (iy = nylo_in; iy <= nyhi_in; iy++)
        for (ix = nxlo_in; ix <= nxhi_in; ix++) {
          variable_electronic_specific_heat =
              cetable_active ? linearinterpolation(T_electron_old[iz][iy][ix], "ce")
                             : electronic_specific_heat;                             
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
            inner_dt/(variable_electronic_specific_heat*electronic_density) *
            (heat_flux_gradient(ix, iy, iz, dxinv, dyinv, dzinv) -
             net_energy_transfer[iz][iy][ix]/volgrid);
        }

    // communicate new T_electron values to ghost grid points

    grid->forward_comm(Grid3d::FIX,this,0,1,sizeof(double),
                       grid_buf1,grid_buf2,MPI_DOUBLE);
  }
}

/* ----------------------------------------------------------------------
   read electron temperatures on grid from a user-specified file
------------------------------------------------------------------------- */

void FixTTMCascade::read_electron_temperatures(const std::string &filename)
{
  memory->create3d_offset(T_electron_read, nzlo_in, nzhi_in, nylo_in, nyhi_in, nxlo_in, nxhi_in,
                          "ttm/cascade:T_electron_read");
  memset(&T_electron_read[nzlo_in][nylo_in][nxlo_in], 0, ngridown * sizeof(int));

  // proc 0 opens file

  SafeFilePtr fp;
  if (comm->me == 0) {
    fp = utils::open_potential(filename, lmp, nullptr);
    if (!fp) error->one(FLERR, "Cannot open grid file: {}: {}", filename,
                 utils::getsyserror());
  }

  // read the file
  // Grid3d::read_file() calls back to unpack_read_grid() with chunks of lines

  grid->read_file(Grid3d::FIX,this,fp,CHUNK,MAXLINE);

  // check completeness of input data

  int flag = 0;
  double sum = 0.0;
  for (int iz = nzlo_in; iz <= nzhi_in; iz++)
    for (int iy = nylo_in; iy <= nyhi_in; iy++)
      for (int ix = nxlo_in; ix <= nxhi_in; ix++) {
        if (T_electron_read[iz][iy][ix] == 0) flag = 1;
        sum += T_electron[iz][iy][ix];
      }

  int flagall;
  MPI_Allreduce(&flag, &flagall, 1, MPI_INT, MPI_SUM, world);
  if (flagall) error->all(FLERR, "Fix ttm/cascade infile did not set all temperatures");

  memory->destroy3d_offset(T_electron_read, nzlo_in, nylo_in, nxlo_in);

  double global_sum;
  MPI_Allreduce(&sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, world);
  average_electronic_temperature = global_sum / ngridtotal;

  if (cetable_active)
    electronic_specific_heat = linearinterpolation(
        average_electronic_temperature, "ce"); // To check stability_criterion

  if (ketable_active)
    electronic_thermal_conductivity = linearinterpolation(
        average_electronic_temperature, "ke"); // To check stability_criterion

  MPI_Bcast(&average_electronic_temperature, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&electronic_specific_heat, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&electronic_thermal_conductivity, 1, MPI_DOUBLE, 0, world);

}

/* ----------------------------------------------------------------------
   process a chunk of lines in buffer
   each proc stores values for grid points it owns
   called back to from Grid3d::read_file()
------------------------------------------------------------------------- */

int FixTTMCascade::unpack_read_grid(int /*nlines*/, char *buffer)
{
  // loop over chunk of lines of grid point values
  // skip comment lines
  // tokenize the line into ix,iy,iz grid index plus temperature value
  // if I own grid point, store the value

  int nread = 0;

  for (const auto &line : utils::split_lines(buffer)) {
    try {
      ValueTokenizer values(utils::trim_comment(line));
      if (values.count() == 0) {
        ;    // ignore comment only or blank lines
      } else if (values.count() == 4) {
        ++nread;

        int ix = values.next_int() - 1;
        int iy = values.next_int() - 1;
        int iz = values.next_int() - 1;

        if (ix < 0 || ix >= nxgrid || iy < 0 || iy >= nygrid || iz < 0 || iz >= nzgrid)
          throw TokenizerException("Fix ttm/cascade invalid grid index in input", "");

        if (ix >= nxlo_in && ix <= nxhi_in && iy >= nylo_in && iy <= nyhi_in && iz >= nzlo_in &&
            iz <= nzhi_in) {
          T_electron[iz][iy][ix] = values.next_double();
          T_electron_read[iz][iy][ix] = 1;
        }
      } else {
        throw TokenizerException("Incorrect format in fix ttm electron grid file", "");
      }
    } catch (std::exception &e) {
      error->one(FLERR, e.what());
    }
  }

  return nread;
}

/* ----------------------------------------------------------------------
   pack state of Fix into one write, but not per-grid values
------------------------------------------------------------------------- */

void FixTTMCascade::write_restart(FILE *fp)
{
  double rlist[4];

  rlist[0] = nxgrid;
  rlist[1] = nygrid;
  rlist[2] = nzgrid;
  rlist[3] = seed;

  if (comm->me == 0) {
    int size = 4 * sizeof(double);
    fwrite(&size,sizeof(int),1,fp);
    fwrite(rlist,sizeof(double),4,fp);
  }
}

/* ----------------------------------------------------------------------
   use state info from restart file to restart the Fix
------------------------------------------------------------------------- */

void FixTTMCascade::restart(char *buf)
{
  auto *rlist = (double *) buf;

  // check that restart grid size is same as current grid size

  int nxgrid_old = static_cast<int> (rlist[0]);
  int nygrid_old = static_cast<int> (rlist[1]);
  int nzgrid_old = static_cast<int> (rlist[2]);

  if (nxgrid_old != nxgrid || nygrid_old != nygrid || nzgrid_old != nzgrid)
    error->all(FLERR,"Must restart fix ttm with same grid size");

  // change RN seed from initial seed, to avoid same Langevin factors
  // just increment by 1, since for RanMars that is a new RN stream

  seed = static_cast<int> (rlist[3]) + 1;
  delete random;
  random = new RanMars(lmp,seed+comm->me);
}

/* ----------------------------------------------------------------------
   write electron temperatures on grid to file
   identical format to infile option, so info can be read in when restarting
   each proc contributes info for its portion of grid
------------------------------------------------------------------------- */

void FixTTMCascade::write_restart_file(const char *file)
{
  // proc 0 opens file and writes header

  if (comm->me == 0) {
    auto outfile = std::string(file) + ".ttm";
    fpout = fopen(outfile.c_str(),"w");
    if (fpout == nullptr)
      error->one(FLERR,"Cannot open fix ttm/cascade restart file {}: {}",outfile,utils::getsyserror());

    utils::print(fpout,"# DATE: {} UNITS: {} COMMENT: "
               "Electron temperature on {}x{}x{} grid at step {} - "
               "created by fix {}\n",
               utils::current_date(),update->unit_style,
               nxgrid,nygrid,nzgrid,update->ntimestep,style);
  }

  // write file
  // Grid3d::write_file() calls back to pack_write_grid() and unpack_write_grid()

  grid->write_file(Grid3d::FIX,this,0,1,sizeof(double), MPI_DOUBLE);

  // close file

  if (comm->me == 0) fclose(fpout);
}

/* ----------------------------------------------------------------------
   pack values from local grid into buf
------------------------------------------------------------------------- */

void FixTTMCascade::pack_write_grid(int /*which*/, void *vbuf)
{
  int ix, iy, iz;

  auto *buf = (double *) vbuf;

  int m = 0;
  for (iz = nzlo_in; iz <= nzhi_in; iz++)
    for (iy = nylo_in; iy <= nyhi_in; iy++)
      for (ix = nxlo_in; ix <= nxhi_in; ix++)
        buf[m++] = T_electron[iz][iy][ix];
}

/* ----------------------------------------------------------------------
   unpack values from buf and write them to restart file
------------------------------------------------------------------------- */

void FixTTMCascade::unpack_write_grid(int /*which*/, void *vbuf, int *bounds)
{
  int ix, iy, iz;

  int xlo = bounds[0];
  int xhi = bounds[1];
  int ylo = bounds[2];
  int yhi = bounds[3];
  int zlo = bounds[4];
  int zhi = bounds[5];

  auto *buf = (double *) vbuf;
  double value;

  int m = 0;
  for (iz = zlo; iz <= zhi; iz++)
    for (iy = ylo; iy <= yhi; iy++)
      for (ix = xlo; ix <= xhi; ix++) {
        value = buf[m++];
        fprintf(fpout, "%d %d %d %20.16g\n", ix+1, iy+1, iz+1, value);
      }
}

/* ----------------------------------------------------------------------
   subset of grid assigned to each proc may have changed
   called by load balancer when proc subdomains are adjusted
------------------------------------------------------------------------- */

void FixTTMCascade::reset_grid()
{
  // check if new grid partitioning is different on any proc
  // if not, just return

  int tmp[12];
  double maxdist = 0.5 * neighbor->skin;
  auto *gridnew = new Grid3d(lmp, world, nxgrid, nygrid, nzgrid);
  gridnew->set_distance(maxdist);
  gridnew->set_stencil_grid(1,1);
  gridnew->setup_grid(tmp[0],tmp[1],tmp[2],tmp[3],tmp[4],tmp[5],
                      tmp[6],tmp[7],tmp[8],tmp[9],tmp[10],tmp[11]);

  if (grid->identical(gridnew)) {
    delete gridnew;
    return;
  } else delete gridnew;

  // delete grid data which doesn't need to persist from previous to new decomp

  memory->destroy(grid_buf1);
  memory->destroy(grid_buf2);
  memory->destroy3d_offset(T_electron_old, nzlo_out, nylo_out, nxlo_out);
  memory->destroy3d_offset(net_energy_transfer, nzlo_out, nylo_out, nxlo_out);

  // make copy of ptrs to grid data which does need to persist

  grid_previous = grid;
  T_electron_previous = T_electron;
  nxlo_out_previous = nxlo_out;
  nylo_out_previous = nylo_out;
  nzlo_out_previous = nzlo_out;

  // allocate new per-grid data for new decomposition

  allocate_grid();

  // perform remap from previous decomp to new decomp

  int nremap_buf1,nremap_buf2;
  grid->setup_remap(grid_previous,nremap_buf1,nremap_buf2);

  double *remap_buf1,*remap_buf2;
  memory->create(remap_buf1, nremap_buf1, "ttm/cascade:remap_buf1");
  memory->create(remap_buf2, nremap_buf2, "ttm/cascade:remap_buf2");

  grid->remap(Grid3d::FIX,this,0,1,sizeof(double),remap_buf1,remap_buf2,MPI_DOUBLE);

  memory->destroy(remap_buf1);
  memory->destroy(remap_buf2);

  // delete grid data and grid for previous decomposition

  memory->destroy3d_offset(T_electron_previous,
                           nzlo_out_previous, nylo_out_previous,
                           nxlo_out_previous);
  delete grid_previous;

  // communicate temperatures to ghost cells on new grid

  grid->forward_comm(Grid3d::FIX,this,0,1,sizeof(double),
                     grid_buf1,grid_buf2,MPI_DOUBLE);

  // zero new net_energy_transfer
  // in case compute_vector accesses it on timestep 0

  outflag = 0;
  memset(&net_energy_transfer[nzlo_out][nylo_out][nxlo_out],0,
         ngridout*sizeof(double));
}

/* ----------------------------------------------------------------------
   pack own values to buf to send to another proc
------------------------------------------------------------------------- */

void FixTTMCascade::pack_forward_grid(int /*which*/, void *vbuf, int nlist, int *list)
{
  auto *buf = (double *) vbuf;
  double *src = &T_electron[nzlo_out][nylo_out][nxlo_out];

  for (int i = 0; i < nlist; i++) buf[i] = src[list[i]];
}

/* ----------------------------------------------------------------------
   unpack another proc's own values from buf and set own ghost values
------------------------------------------------------------------------- */

void FixTTMCascade::unpack_forward_grid(int /*which*/, void *vbuf, int nlist, int *list)
{
  auto *buf = (double *) vbuf;
  double *dest = &T_electron[nzlo_out][nylo_out][nxlo_out];

  for (int i = 0; i < nlist; i++) dest[list[i]] = buf[i];
}

/* ----------------------------------------------------------------------
   pack ghost values into buf to send to another proc
------------------------------------------------------------------------- */

void FixTTMCascade::pack_reverse_grid(int /*which*/, void *vbuf, int nlist, int *list)
{
  auto *buf = (double *) vbuf;
  double *src = &net_energy_transfer[nzlo_out][nylo_out][nxlo_out];

  for (int i = 0; i < nlist; i++) buf[i] = src[list[i]];
}

/* ----------------------------------------------------------------------
   unpack another proc's ghost values from buf and add to own values
------------------------------------------------------------------------- */

void FixTTMCascade::unpack_reverse_grid(int /*which*/, void *vbuf, int nlist, int *list)
{
  auto *buf = (double *) vbuf;
  double *dest = &net_energy_transfer[nzlo_out][nylo_out][nxlo_out];

  for (int i = 0; i < nlist; i++) dest[list[i]] += buf[i];
}

/* ----------------------------------------------------------------------
   pack old grid values to buf to send to another proc
------------------------------------------------------------------------- */

void FixTTMCascade::pack_remap_grid(int /*which*/, void *vbuf, int nlist, int *list)
{
  auto *buf = (double *) vbuf;
  double *src =
    &T_electron_previous[nzlo_out_previous][nylo_out_previous][nxlo_out_previous];

  for (int i = 0; i < nlist; i++) buf[i] = src[list[i]];
}

/* ----------------------------------------------------------------------
   unpack another proc's own values from buf and set own ghost values
------------------------------------------------------------------------- */

void FixTTMCascade::unpack_remap_grid(int /*which*/, void *vbuf, int nlist, int *list)
{
  auto *buf = (double *) vbuf;
  double *dest = &T_electron[nzlo_out][nylo_out][nxlo_out];

  for (int i = 0; i < nlist; i++) dest[list[i]] = buf[i];
}

/* ----------------------------------------------------------------------
   allocate 3d grid quantities
------------------------------------------------------------------------- */

void FixTTMCascade::allocate_grid()
{
  double maxdist = 0.5 * neighbor->skin;

  grid = new Grid3d(lmp, world, nxgrid, nygrid, nzgrid);
  grid->set_distance(maxdist);
  grid->set_stencil_grid(1,1);
  grid->setup_grid(nxlo_in, nxhi_in, nylo_in, nyhi_in, nzlo_in, nzhi_in,
                   nxlo_out, nxhi_out, nylo_out, nyhi_out, nzlo_out, nzhi_out);

  ngridown = (nxhi_in - nxlo_in + 1) * (nyhi_in - nylo_in + 1) *
    (nzhi_in - nzlo_in + 1);
  ngridout = (nxhi_out - nxlo_out + 1) * (nyhi_out - nylo_out + 1) *
    (nzhi_out - nzlo_out + 1);

  // setup grid communication and allocate grid data structs

  grid->setup_comm(ngrid_buf1, ngrid_buf2);

  memory->create(grid_buf1, ngrid_buf1, "ttm/cascade:grid_buf1");
  memory->create(grid_buf2, ngrid_buf2, "ttm/cascade:grid_buf2");

  memory->create3d_offset(T_electron_old, nzlo_out, nzhi_out, nylo_out, nyhi_out, nxlo_out,
                          nxhi_out, "ttm/cascade:T_electron_old");
  memory->create3d_offset(T_electron, nzlo_out, nzhi_out, nylo_out, nyhi_out, nxlo_out, nxhi_out,
                          "ttm/cascade:T_electron");
  memory->create3d_offset(net_energy_transfer, nzlo_out, nzhi_out, nylo_out, nyhi_out, nxlo_out,
                          nxhi_out, "ttm/cascade:net_energy_transfer");
  memory->create3d_offset(conductivity_xface, nzlo_out, nzhi_out, nylo_out, nyhi_out, nxlo_out,
                          nxhi_out, "ttm/cascade:conductivity_xface");
  memory->create3d_offset(conductivity_yface, nzlo_out, nzhi_out, nylo_out, nyhi_out, nxlo_out,
                          nxhi_out, "ttm/cascade:conductivity_yface");
  memory->create3d_offset(conductivity_zface, nzlo_out, nzhi_out, nylo_out, nyhi_out, nxlo_out,
                          nxhi_out, "ttm/cascade:conductivity_zface");
}

/* ----------------------------------------------------------------------
   deallocate 3d grid quantities
------------------------------------------------------------------------- */

void FixTTMCascade::deallocate_grid()
{
  delete grid;
  memory->destroy(grid_buf1);
  memory->destroy(grid_buf2);

  memory->destroy3d_offset(T_electron_old, nzlo_out, nylo_out, nxlo_out);
  memory->destroy3d_offset(T_electron, nzlo_out, nylo_out, nxlo_out);
  memory->destroy3d_offset(net_energy_transfer, nzlo_out, nylo_out, nxlo_out);
  memory->destroy3d_offset(conductivity_xface, nzlo_out, nylo_out, nxlo_out);
  memory->destroy3d_offset(conductivity_yface, nzlo_out, nylo_out, nxlo_out);
  memory->destroy3d_offset(conductivity_zface, nzlo_out, nylo_out, nxlo_out);

}

/* ----------------------------------------------------------------------
   return index of grid associated with name
   this class can store M named grids, indexed 0 to M-1
   also set dim for 2d vs 3d grid
   return -1 if grid name not found
------------------------------------------------------------------------- */

int FixTTMCascade::get_grid_by_name(const std::string &name, int &dim)
{
  if (name == "grid") {
    dim = 3;
    return 0;
  }

  return -1;
}

/* ----------------------------------------------------------------------
   return ptr to Grid data struct for grid with index
   this class can store M named grids, indexed 0 to M-1
   return nullptr if index is invalid
------------------------------------------------------------------------- */

void *FixTTMCascade::get_grid_by_index(int index)
{
  if (index == 0) return grid;

  return nullptr;
}

/* ----------------------------------------------------------------------
   return index of data associated with name in grid with index igrid
   this class can store M named grids, indexed 0 to M-1
   each grid can store G named data sets, indexed 0 to G-1
     a data set name can be associated with multiple grids
   set ncol for data set, 0 = vector, 1-N for array with N columns
     vector = single value per grid pt, array = N values per grid pt
   return -1 if data name not found
------------------------------------------------------------------------- */

int FixTTMCascade::get_griddata_by_name(int igrid, const std::string &name, int &ncol)
{
  if ((igrid == 0) && (name == "data")) {
    ncol = 0;
    return 0;
  }

  return -1;
}

/* ----------------------------------------------------------------------
   return ptr to multidim data array associated with index
   this class can store G named data sets, indexed 0 to M-1
   return nullptr if index is invalid
------------------------------------------------------------------------- */

void *FixTTMCascade::get_griddata_by_index(int index)
{
  if (index == 0) return T_electron;

  return nullptr;
}

/* ----------------------------------------------------------------------
   return the energy of the electronic subsystem
   or the net_energy transfer between the subsystems
------------------------------------------------------------------------- */

double FixTTMCascade::compute_vector(int n) {
  int ix, iy, iz;

  if (outflag == 0) {
    double dx = domain->xprd / nxgrid;
    double dy = domain->yprd / nygrid;
    double dz = domain->zprd / nzgrid;
    double volgrid = dx * dy * dz;

    double e_energy_me = 0.0;
    double transfer_energy_me = 0.0;

    for (iz = nzlo_in; iz <= nzhi_in; iz++)
      for (iy = nylo_in; iy <= nyhi_in; iy++)
        for (ix = nxlo_in; ix <= nxhi_in; ix++) {

          e_energy_me += integrated_ce(T_electron[iz][iy][ix]) *
                         electronic_density * volgrid;
          transfer_energy_me += net_energy_transfer[iz][iy][ix] * update->dt;
        }

    MPI_Allreduce(&e_energy_me, &e_energy, 1, MPI_DOUBLE, MPI_SUM, world);
    MPI_Allreduce(&transfer_energy_me, &transfer_energy, 1, MPI_DOUBLE, MPI_SUM, world);
    outflag = 1;
  }

  if (n == 0) return e_energy;
  if (n == 1) return transfer_energy;
  return 0.0;
}

/* ----------------------------------------------------------------------
   memory usage for flangevin and 3d grid
------------------------------------------------------------------------- */

double FixTTMCascade::memory_usage()
{
  double bytes = 0.0;
  bytes += (double) 3 * atom->nmax * sizeof(double);
  bytes += (double) 3 * ngridout * sizeof(double);
  return bytes;
}

/* ----------------------------------------------------------------------
   read the specific heat and thermal conductivity table files
------------------------------------------------------------------------- */

void FixTTMCascade::tableinterpreader(const std::string &filename,
                                   const std::string &keyword) {

  std::vector<double> &temp_vals =
      (keyword == "ce") ? temp_ce_values : temp_ke_values;
  std::vector<double> &dtemp_vals =
      (keyword == "ce") ? dtemp_ce_values : dtemp_ke_values;
  std::vector<double> &y_vals = (keyword == "ce") ? ce_values : ke_values;
  std::vector<double> &dy_vals = (keyword == "ce") ? dce_values : dke_values;

  bool table_flag = false;

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
    if (keyword == "ce") {
      ce_integral_values.resize(nsize_table);
      ce_integral_values[0] = y_vals[0] * temp_vals[0];
    }
    for (int i = 0; i < nsize_table - 1; i++) {
      dtemp_vals[i] = temp_vals[i + 1] - temp_vals[i];
      dy_vals[i] = y_vals[i + 1] - y_vals[i];
      if (keyword == "ce") {
        // numerical integration
        ce_integral_values[i + 1] = ce_integral_values[i] + 0.5 * (y_vals[i + 1] + y_vals[i]) * (temp_vals[i + 1] - temp_vals[i]);
      }
      if (temp_vals[i] >= temp_vals[i + 1]) table_flag = true;
    }
  }


  int nsize_table = static_cast<int>(temp_vals.size());
  MPI_Bcast(&nsize_table, 1, MPI_INT, 0, world);

  if (comm->me != 0) {
    temp_vals.resize(nsize_table);
    y_vals.resize(nsize_table);
    dtemp_vals.resize(nsize_table);
    dy_vals.resize(nsize_table);
    if (keyword == "ce")
      ce_integral_values.resize(nsize_table);
  }

  MPI_Bcast(temp_vals.data(), nsize_table, MPI_DOUBLE, 0, world);
  MPI_Bcast(y_vals.data(), nsize_table, MPI_DOUBLE, 0, world);
  MPI_Bcast(dtemp_vals.data(), nsize_table, MPI_DOUBLE, 0, world);
  MPI_Bcast(dy_vals.data(), nsize_table, MPI_DOUBLE, 0, world);
  MPI_Bcast(&table_flag, 1, MPI_C_BOOL, 0, world);
  if (keyword == "ce") MPI_Bcast(ce_integral_values.data(), nsize_table, MPI_DOUBLE, 0, world);

  // error check

  if (table_flag) error->all(FLERR, "Two consecutive values in the {} table are the same or not sorted", keyword);

  if (temp_vals.size() < 2) error->all(FLERR, "The {} table has less than 2 rows, interpolation is invalid", keyword);


  }

/* ----------------------------------------------------------------------
   interpolates the specific heat and thermal conductivity values
------------------------------------------------------------------------- */

double FixTTMCascade::linearinterpolation(double temp, const std::string &keyword) {

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

/* ----------------------------------------------------------------------
   performs the numerical integration of the specific heat
------------------------------------------------------------------------- */

double FixTTMCascade::integrated_ce(double Te) {
  if (!cetable_active) {
    return electronic_specific_heat * Te;
  }

  int lo = 0;
  int hi = static_cast<int>(ce_values.size()) - 1;
  int mid;

  if (Te <= temp_ce_values[lo])
    // Exact match or lower than lower bound
    return ce_values[lo] * Te;
  if (Te >= temp_ce_values[hi])
    // Exact match or higher than upper bound
    return ce_integral_values[hi] + (ce_values[hi] * (Te - temp_ce_values[hi]));

  // Look for fraction
  while (hi - lo > 1) {
    mid = (lo + hi) / 2;
    if (Te < temp_ce_values[mid])
      hi = mid;
    else
      lo = mid;
  }

  double Ce_interpolated =
      ce_values[lo] +
      (Te - temp_ce_values[lo]) * (dce_values[lo] / dtemp_ce_values[lo]);

  return ce_integral_values[lo] +
         0.5 * (ce_values[lo] + Ce_interpolated) * (Te - temp_ce_values[lo]);
}

/* ----------------------------------------------------------------------
   computes the heat flux gradient according to the type of thermal conductivity.
   constant ke: k \nabla^2 T
   variable ke: \nabla * q, where q = k * \nabla T
------------------------------------------------------------------------- */

double FixTTMCascade::heat_flux_gradient(int ix, int iy, int iz, double dxinv, double dyinv,
                       double dzinv) {
  
  double heat_flux_gradient;

  if (!ketable_active) {

    heat_flux_gradient = electronic_thermal_conductivity * ((T_electron_old[iz][iy][ix-1] + T_electron_old[iz][iy][ix+1] -
               2.0*T_electron_old[iz][iy][ix])*dxinv*dxinv +
              (T_electron_old[iz][iy-1][ix] + T_electron_old[iz][iy+1][ix] -
               2.0*T_electron_old[iz][iy][ix])*dyinv*dyinv +
              (T_electron_old[iz-1][iy][ix] + T_electron_old[iz+1][iy][ix] -
               2.0*T_electron_old[iz][iy][ix])*dzinv*dzinv);

    return heat_flux_gradient;           
  }                      

  else {

    heat_flux_gradient = 
              (conductivity_xface[iz][iy][ix] * (T_electron_old[iz][iy][ix+1] - T_electron_old[iz][iy][ix]) // dq/dx
              + conductivity_xface[iz][iy][ix-1] * (T_electron_old[iz][iy][ix-1] - T_electron_old[iz][iy][ix]))*dxinv*dxinv
              + (conductivity_yface[iz][iy][ix] * (T_electron_old[iz][iy+1][ix] - T_electron_old[iz][iy][ix]) // dq/dy
              + conductivity_yface[iz][iy-1][ix] * (T_electron_old[iz][iy-1][ix] - T_electron_old[iz][iy][ix]))*dyinv*dyinv
              + (conductivity_zface[iz][iy][ix] * (T_electron_old[iz+1][iy][ix] - T_electron_old[iz][iy][ix]) // dq/dz
              + conductivity_zface[iz-1][iy][ix] * (T_electron_old[iz-1][iy][ix] - T_electron_old[iz][iy][ix]))*dzinv*dzinv;

    return heat_flux_gradient;            
  }  
}
