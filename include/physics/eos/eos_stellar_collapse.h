/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2018 Triad National Security, LLC
 * All rights reserved.
 * --------------------------------------------------------------------------~*/

/******************************************************************************
 *                                                                            *
 * EOS_STELLAR_COLLAPSE.c                                                     *
 * TODO: CLEANUP                                                              *
 *                                                                            *
 * IMPLEMENTS ROUTINES FOR READING EOS TABLES PROVIDED ON STELLARCOLLAPSE.ORG *
 *                                                                            *
 ******************************************************************************/
#pragma once

#include "params.h"
#include <fstream>

#if 0
// TODO : Need to coporate with FleCSPH unit system
//        below const definition doesn't change unit.
//        just define like that to work it
const double RHO_unit = 1.0; // For density
const double U_unit = 1.0;  // For internel specific energy
#endif

// HDF5
#include <hdf5.h>

using std::isnan;
using std::isinf;

namespace eos {

template<>
class eos_t<param::eos_stellar_collapse>{
public:
  /**
  * @brief      Initialize tabulated EOS from stellarcollapse
  *             Uses the path to EOS table (in HDF5 format).
  */
  static void init() {
    log_one(info) << "Reading tabulated EOS from file: "
                << param::eos_tab_file_path << std::endl;
    EOS_SC_init(param::eos_tab_file_path);
  }


  /**
  * @brief      Compute pressure for tabulated EOS
  * @param      particle
  */
  static void compute_pressure(body & particle) {
    double pressure = EOS_pressure_rho0_u(particle);
    particle.setPressure(pressure);
  } // compute_pressure_sc


  /**
  * @brief      Compute speed of sound for tabulated EOS
  * @param      particle
  */
  static void
  compute_soundspeed(body & particle) {
    double soundspeed = EOS_sound_speed_rho0_u(particle);
    particle.setSoundspeed(soundspeed);
  } // compute_soundspeed_sc

  /**
  * @brief      Compute entropy
  *             TODO: implement
  *
  * @param      particle
  */
  static void
  compute_entropy(body & particle) {
    /* ... */
  }

  /**
  * @brief      Compute temperature for tabulated EOS
  * @param      particle
  */
  static void
  compute_temperature(body & particle) {
    double temperature = EOS_temperature_sc(particle);
    particle.setTemperature(temperature);
  } // compute_temperature_sc

  static void compute_internal_energy(body& particle){
    const double MEV = 1.60217653e-6, // [erg/MeV] - conversion factor
      KBOL = 1.3806505e-16; // [erg/K]
    const double rho = particle.getDensity(),
               T = particle.getTemperature() * KBOL / MEV, // T in MeV
    ye = particle.getElectronfraction();
    double u = eos_t<param::eos_stellar_collapse>::EOS_SC_get_u_of_T(rho, T, ye);
    particle.setInternalenergy(u);
  }

  /**
  * @brief      Compute entropy, pressure, soundspeed and temperature
  *             TODO: needs more work
  *
  * @param      particle
  */
  static void
  compute_spct_given_rho_u(body & particle) {
    compute_entropy(particle);
    compute_pressure(particle);
    compute_soundspeed(particle);
    compute_temperature(particle);
  }

  // TODO
  static get_quantity_t get_dpdrho_at_temp;

private:

  // fail = 0 = secant or bisection failed
  // degenerate = 1 = root find succeeded, but wrong root
  // success_boundary = 2 = boundaries of the table for plotting
  enum record_type { fail, degen, success_boundary };

  // struct to store a particular failed root-find
  struct table_val {
    double lrho; // log density
    double lT; // log temperature
    double Ye; // log electron per baryon fraction
    double le; // log specific internal energy -- shifted
    double lT_root; // root find attempt
    record_type rt;
    int lr_idx;
    int Ye_idx;
    // std::vector<double> lTs; // All temperatures for this lrho,Ye point
    // std::vector<double> les; // All energies = f(lT) for this lrho,Ye point
  };

  struct of_lT_params {
    double lrho, ye;
    double * tab;
  };
  struct of_lT_adiabat_params {
    double * tab;
    const struct of_adiabat * a;
  };

  static constexpr double TABLE_TOL = 1.e-10;
  static constexpr double TABLE_FTOL = 1.e-10;
  static constexpr int SC_DEBUG = 0;
  static constexpr int SC_MONOTONE_SAFE = 1;
  static constexpr int SC_THROTTLE_CS = 0;

  static auto
  EOS_ELEM(int irho, int iT, int iY) {
    return (Nrho * ((iY)*NT + (iT)) + (irho));
  }
  static auto
  YE_ELEM(int i, int j, int k) {
    return (NYe_ye * ((i)*NT_ye + (j)) + (k));
  }

  static auto
  MMA_ELEM(int irho, int iY){
    return (Nrho * iY + irho);
  }

  static int Nrho, NT, NYe, Nrho_ye, NT_ye, NYe_ye;
  static double * tab_lrho;
  static double * tab_lT;
  static double * tab_Ye;
  static double * tab_lP;
  static double * tab_ent;
  static double * tab_dpderho;
  static double * tab_dpdrhoe;
  static double * tab_cs2;
  static double * tab_le;
  static double * tab_Xa;
  static double * tab_Xh;
  static double * tab_Xn;
  static double * tab_Xp;
  static double * tab_Abar;
  static double * tab_Zbar;
  static double * tab_lwmrho; // log enthalpy - rho, by volume
  static double * tab_hm1; // enthalpy - 1, by mass
  static double * tab_poly_gamma; // Polytrope gamma
  static double * tab_poly_K; // polytrope K
  static double * tab_rho;
  static double * tab_Yeye;
  static double * tab_T;
  static double * tab_dYedt;
  static double * tab_deweakdt;

  // min and max of wmrho given fixed ilrho and iY
  static double * tab_le_min_2d;
  static double * tab_le_max_2d;
  static double * tab_lP_min_2d;
  static double * tab_lP_max_2d;
  static double * tab_lwmrho_min_2d;
  static double * tab_lwmrho_max_2d;
  static double * tab_hm1_min_1d;

  static double tab_lrho_min, tab_lrho_max;
  static double tab_rhoye_min, tab_rhoye_max;
  static double tab_lT_min, tab_lT_max;
  static double tab_Tye_min, tab_Tye_max;
  static double tab_Ye_min, tab_Ye_max;
  static double tab_Yeye_min, tab_Yeye_max;
  static double tab_dlrho, tab_dlT, tab_dYe;
  static double tab_drhoye, tab_dTye, tab_dYeye;

  static double tab_lP_min, tab_lP_max;
  static double tab_ent_min, tab_ent_max;
  static double tab_cs2_min, tab_cs2_max;
  static double tab_le_min, tab_le_max;
  static double tab_Xa_min, tab_Xa_max;
  static double tab_Xh_min, tab_Xh_max;
  static double tab_Xn_min, tab_Xn_max;
  static double tab_Xp_min, tab_Xp_max;
  static double tab_Abar_min, tab_Abar_max;
  static double tab_Zbar_min, tab_Zbar_max;
  static double tab_dpderho_min, tab_dpderho_max;
  static double tab_dpdrhoe_min, tab_dpdrhoe_max;
  static double tab_lwmrho_min, tab_lwmrho_max;
  static double tab_dYedt_min, tab_dYedt_max;
  static double tab_deweakdt_min, tab_deweakdt_max;

  static double tab_rho_min, tab_rho_max;
  static double tab_T_min, tab_T_max;
  static double tab_e_min, tab_e_max;
  static double tab_P_min, tab_P_max;
  static double tab_wmrho_min, tab_wmrho_max;
  static double tab_hm1_min, tab_hm1_max;

  static double pressure_min;

  static double energy_shift;
  static double enthalpy_shift;

  // Init
  // ----------------------------------------------------------------------
  static void EOS_SC_init(const char * name) {

    std::ifstream infile(name);
    if(!infile.good()) {
      log_one(error) << "File " << name << " not found." << std::endl;
      MPI_Finalize();
      exit(-1);
    }

    hid_t file_id =
      H5Fopen(name, H5F_ACC_RDONLY, H5P_DEFAULT); // opens existing hdf5 file

    hid_t dset_id = H5Dopen(file_id, "pointsrho", H5P_DEFAULT);
    int status =
      H5Dread(dset_id, H5T_STD_I32LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &Nrho);
    H5Dclose(dset_id);

    dset_id = H5Dopen(file_id, "pointstemp", H5P_DEFAULT);
    status = H5Dread(dset_id, H5T_STD_I32LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &NT);
    H5Dclose(dset_id);

    dset_id = H5Dopen(file_id, "pointsye", H5P_DEFAULT);
    status = H5Dread(dset_id, H5T_STD_I32LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &NYe);
    H5Dclose(dset_id);

    dset_id = H5Dopen(file_id, "energy_shift", H5P_DEFAULT);
    status = H5Dread(
      dset_id, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &energy_shift);
    H5Dclose(dset_id);

    int tab_size = NYe * NT * Nrho;
    int tab_size_2d = NYe * Nrho;

    tab_lrho = safe_malloc<double>(Nrho);
    tab_lT = safe_malloc<double>(NT);
    tab_Ye = safe_malloc<double>(NYe);

    tab_lP = safe_malloc<double>(tab_size);
    tab_ent = safe_malloc<double>(tab_size);
    tab_cs2 = safe_malloc<double>(tab_size);
    tab_Xa = safe_malloc<double>(tab_size);
    tab_Xh = safe_malloc<double>(tab_size);
    tab_Xn = safe_malloc<double>(tab_size);
    tab_Xp = safe_malloc<double>(tab_size);
    tab_Abar = safe_malloc<double>(tab_size);
    tab_Zbar = safe_malloc<double>(tab_size);
    tab_le = safe_malloc<double>(tab_size);
    tab_lwmrho = safe_malloc<double>(tab_size);
    tab_hm1 = safe_malloc<double>(tab_size);
    tab_dpderho = safe_malloc<double>(tab_size);
    tab_dpdrhoe = safe_malloc<double>(tab_size);
    tab_poly_gamma = safe_malloc<double>(tab_size);
    tab_poly_K = safe_malloc<double>(tab_size);

    tab_le_min_2d = safe_malloc<double>(tab_size_2d);
    tab_le_max_2d = safe_malloc<double>(tab_size_2d);
    tab_lP_min_2d = safe_malloc<double>(tab_size_2d);
    tab_lP_max_2d = safe_malloc<double>(tab_size_2d);
    tab_lwmrho_min_2d = safe_malloc<double>(tab_size_2d);
    tab_lwmrho_max_2d = safe_malloc<double>(tab_size_2d);
    tab_hm1_min_1d = safe_malloc<double>(NYe);

    dset_id = H5Dopen(file_id, "logrho", H5P_DEFAULT);
    status =
      H5Dread(dset_id, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, tab_lrho);
    H5Dclose(dset_id);

    dset_id = H5Dopen(file_id, "logtemp", H5P_DEFAULT);
    status =
      H5Dread(dset_id, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, tab_lT);
    H5Dclose(dset_id);

    dset_id = H5Dopen(file_id, "ye", H5P_DEFAULT);
    status =
      H5Dread(dset_id, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, tab_Ye);
    H5Dclose(dset_id);

    dset_id = H5Dopen(file_id, "logpress", H5P_DEFAULT);
    status =
      H5Dread(dset_id, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, tab_lP);
    H5Dclose(dset_id);

    dset_id = H5Dopen(file_id, "logenergy", H5P_DEFAULT);
    status =
      H5Dread(dset_id, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, tab_le);
    H5Dclose(dset_id);

    dset_id = H5Dopen(file_id, "entropy", H5P_DEFAULT);
    status =
      H5Dread(dset_id, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, tab_ent);
    H5Dclose(dset_id);

    dset_id = H5Dopen(file_id, "dpdrhoe", H5P_DEFAULT);
    status = H5Dread(
      dset_id, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, tab_dpdrhoe);
    H5Dclose(dset_id);

    dset_id = H5Dopen(file_id, "dpderho", H5P_DEFAULT);
    status = H5Dread(
      dset_id, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, tab_dpderho);
    H5Dclose(dset_id);

    dset_id = H5Dopen(file_id, "Xa", H5P_DEFAULT);
    status =
      H5Dread(dset_id, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, tab_Xa);
    H5Dclose(dset_id);

    dset_id = H5Dopen(file_id, "Xh", H5P_DEFAULT);
    status =
      H5Dread(dset_id, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, tab_Xh);
    H5Dclose(dset_id);

    dset_id = H5Dopen(file_id, "Xn", H5P_DEFAULT);
    status =
      H5Dread(dset_id, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, tab_Xn);
    H5Dclose(dset_id);

    dset_id = H5Dopen(file_id, "Xp", H5P_DEFAULT);
    status =
      H5Dread(dset_id, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, tab_Xp);
    H5Dclose(dset_id);

    dset_id = H5Dopen(file_id, "Abar", H5P_DEFAULT);
    status =
      H5Dread(dset_id, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, tab_Abar);
    H5Dclose(dset_id);

    dset_id = H5Dopen(file_id, "Zbar", H5P_DEFAULT);
    status =
      H5Dread(dset_id, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, tab_Zbar);
    H5Dclose(dset_id);

    dset_id = H5Dopen(file_id, "gamma", H5P_DEFAULT);
    status = H5Dread(
      dset_id, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, tab_poly_gamma);
    H5Dclose(dset_id);

    dset_id = H5Dopen(file_id, "cs2", H5P_DEFAULT);
    status =
      H5Dread(dset_id, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, tab_cs2);
    H5Dclose(dset_id);

    H5Fclose(file_id);

    // mins and maxes of indep. variables
    tab_lrho_min = tab_lrho[0];
    tab_lrho_max = tab_lrho[Nrho - 1];
    tab_dlrho = (tab_lrho_max - tab_lrho_min) / (Nrho - 1);
    tab_rho_min = pow(10., tab_lrho_min);
    tab_rho_max = pow(10., tab_lrho_max);

    tab_lT_min = tab_lT[0];
    tab_lT_max = tab_lT[NT - 1];
    tab_dlT = (tab_lT_max - tab_lT_min) / (NT - 1);
    tab_T_min = pow(10., tab_lT_min);
    tab_T_max = pow(10., tab_lT_max);

    tab_Ye_min = tab_Ye[0];
    tab_Ye_max = tab_Ye[NYe - 1];
    tab_dYe = (tab_Ye_max - tab_Ye_min) / (NYe - 1);

    // mins and maxes of dependent variables
    tab_lP_min = find_min(tab_lP, tab_size);
    tab_lP_max = find_max(tab_lP, tab_size);
    tab_ent_min = find_min(tab_ent, tab_size);
    tab_ent_max = find_max(tab_ent, tab_size);
    tab_le_min = find_min(tab_le, tab_size);
    tab_le_max = find_max(tab_le, tab_size);
    tab_Xa_min = find_min(tab_Xa, tab_size);
    tab_Xa_max = find_max(tab_Xa, tab_size);
    tab_Xh_min = find_min(tab_Xh, tab_size);
    tab_Xh_max = find_max(tab_Xh, tab_size);
    tab_Xn_min = find_min(tab_Xn, tab_size);
    tab_Xn_max = find_max(tab_Xn, tab_size);
    tab_Xp_min = find_min(tab_Xp, tab_size);
    tab_Xp_max = find_max(tab_Xp, tab_size);
    tab_Abar_min = find_min(tab_Abar, tab_size);
    tab_Abar_max = find_max(tab_Abar, tab_size);
    tab_Zbar_min = find_min(tab_Zbar, tab_size);
    tab_Zbar_max = find_max(tab_Zbar, tab_size);
    tab_dpdrhoe_min = find_min(tab_dpdrhoe, tab_size);
    tab_dpdrhoe_max = find_max(tab_dpdrhoe, tab_size);
    tab_dpderho_min = find_min(tab_dpderho, tab_size);
    tab_dpderho_max = find_max(tab_dpderho, tab_size);
    tab_cs2_min = find_min(tab_cs2, tab_size);
    tab_cs2_max = find_max(tab_cs2, tab_size);

    tab_e_min = le2e(tab_le_min);
    tab_e_max = le2e(tab_le_max);
    tab_P_min = pow(10., tab_lP_min);
    tab_P_max = pow(10., tab_lP_max);
    fill_max_min_2d(tab_le_max_2d, tab_le_min_2d, tab_le);
    fill_max_min_2d(tab_lP_max_2d, tab_lP_min_2d, tab_lP);

    { // make enthalpy table
      double * tab_wmrho = safe_malloc<double>(tab_size);
      // double *tab_hm1   = safe_malloc(tab_size*sizeof(double));
      double lrho, rho, le, lP, e, P, w, lw, h;
      // first tabulate the real enthalpy
      for(int irho = 0; irho < Nrho; irho++) {
        lrho = tab_lrho[irho];
        rho = pow(10., lrho);
        for(int iT = 0; iT < NT; iT++) {
          for(int iY = 0; iY < NYe; iY++) {
            le = tab_le[EOS_ELEM(irho, iT, iY)];
            lP = tab_lP[EOS_ELEM(irho, iT, iY)];
            e = le2e(le);
            P = pow(10., lP);
            // use log enthalpy minus rho.
            w = rho * e + P;
            h = e + P / rho;
            tab_wmrho[EOS_ELEM(irho, iT, iY)] = w;
            tab_hm1[EOS_ELEM(irho, iT, iY)] = h;
          }
        }
      }

      // Calculate min and max of wmrho
      tab_wmrho_min = find_min(tab_wmrho, tab_size);
      tab_wmrho_max = find_max(tab_wmrho, tab_size);
      tab_hm1_min = find_min(tab_hm1, tab_size);
      tab_hm1_max = find_max(tab_hm1, tab_size);
      fill_min_1d(tab_hm1_min_1d, tab_hm1);

      // Next check to see if we need a min value
      if(tab_wmrho_min <=
        0) { // arbitrary. Just need to make enthalpy + shift > 0
        enthalpy_shift = -1.01 * tab_wmrho_min;
      }
      else {
        enthalpy_shift = 0.0;
      }

      // Set the log enthalpy
      for(int irho = 0; irho < Nrho; irho++) {
        for(int iT = 0; iT < NT; iT++) {
          for(int iY = 0; iY < NYe; iY++) {
            w = tab_wmrho[EOS_ELEM(irho, iT, iY)];
            h = tab_hm1[EOS_ELEM(irho, iT, iY)];
            lw = w2lw(w);
            tab_lwmrho[EOS_ELEM(irho, iT, iY)] = lw;
          }
        }
      }
      free(tab_wmrho);
      // free(tab_hm1);
    }
    tab_lwmrho_min = find_min(tab_lwmrho, tab_size);
    tab_lwmrho_max = find_max(tab_lwmrho, tab_size);
    fill_max_min_2d(tab_lwmrho_max_2d, tab_lwmrho_min_2d, tab_lwmrho);

    /*{ // make cs2 table
      int elem;
      for (int irho = 0; irho < Nrho; irho++) {
        double lrho = tab_lrho[irho];
        double rho = pow(10.,lrho);
        for (int iT = 0; iT < NT; iT++) {
          for (int iY = 0; iY < NYe; iY++) {
            elem = EOS_ELEM(irho,iT,iY);
            double hm1 = tab_hm1[elem];
            double h = hm1 + C_LIGHT_CGS*C_LIGHT_CGS;
            double lP = tab_lP[elem];
            double P = pow(10.,lP);
            double dpdrhoe = tab_dpdrhoe[elem];
            double dpderho = tab_dpderho[elem];
            double cs2 = (dpdrhoe + (P/(rho*rho))*dpderho)/h;
            tab_cs2[elem] = fabs(cs2);
          }
        }
      }
    }
    tab_cs2_min  = find_min(tab_cs2,  tab_size);
    tab_cs2_max  = find_max(tab_cs2,  tab_size);
    */

    { // make polytrope table
      int elem;
      for(int irho = 0; irho < Nrho; irho++) {
        double lrho = tab_lrho[irho];
        for(int iT = 0; iT < NT; iT++) {
          for(int iY = 0; iY < NYe; iY++) {
            elem = EOS_ELEM(irho, iT, iY);
            double lP = tab_lP[elem];
            double Gam = tab_poly_gamma[elem];
            double lK = lP - Gam * lrho;
            double K = pow(10., lK);
            tab_poly_K[elem] = K;
          }
        }
      }
    }

    // sanity checks
    if(isnan(tab_lwmrho_min)) {
      fprintf(stderr, "[EOS_SC_init]: log enthalpy is nan.\n");
      exit(1);
    }
    if(status) {
      fprintf(
        stderr, "[EOS_SC_init]: HDF5 Returned an error. Status = %d.\n", status);
      exit(1);
    }
  }

  // HL : Disalbe this now
  #if 0
    void do_ye_fixup(int i, int j, int k,
        double pv[NVAR], double pv_prefloor[NVAR])
    {
      pv[YE] = catch_ye(pv[YE]);
    }
  #endif
  // ----------------------------------------------------------------------

  // Front-facing API
  // ----------------------------------------------------------------------
  // void EOS_SC_fill(double* rhoIn, double* uIn, double* yeIn, double* eos)
  static void EOS_SC_fill(body & b, double * eos) {
    double lTguess, leosTemp;
    double lrho, e, le;
    double u = b.getInternalenergy();
    double rho = b.getDensity();
    double ye = b.getElectronfraction();

    // double u      = uIn[UU];
    // double rho    = rhoIn[RHO];
    // double ye     = yeIn[YE];
    // double yedens = p[YE];
    double lT = eos[EOS_LT];
    if(isnan(lT) == true)
      std::cout << "Particle id after lT :      " << b.id() << std::endl;
    // double ye     = yedens / (fabs(rho) + SMALL);

    // into CGS
    u *= GV::U_unit;
    rho *= GV::RHO_unit;

    // dont' fall off the table
    lrho = catch_rho(rho);
    ye = catch_ye(ye);
    e = u;
    le = catch_e(e);
    // le = catch_var_2d(lrho,ye,le,tab_le_min_2d,tab_le_max_2d);

  #if 0
  if constexpr (SC_MONOTONE_SAFE){
      le = catch_var_2d_monotone(lrho,ye,le,tab_le);
  }
  #endif
    // Get a good guess
    lTguess = catch_lT(lT);

  #if 0 // HL : We may need this but not now
      // crash and die if something went wrong here
  if constexpr (SC_DEBUG){
      if (isnan(le)) {
        fprintf(stderr,"[EOS_SC_fill %d]: NAN detected!\n",mpi_io_proc());
        fprintf(stderr,"rho     = %.10e\n",rho);
        fprintf(stderr,"u       = %.10e\n",u);
        fprintf(stderr,"e       = %.10e\n",e);
        fprintf(stderr,"lrho    = %.10f\n",lrho);
        fprintf(stderr,"ye      = %.10f\n",ye);
        fprintf(stderr,"le      = %.10f\n",le);
        fprintf(stderr,"lTguess = %.10f\n",lTguess);
        exit(1);
      }
  }
  #endif

    int status = find_lT(lrho, lTguess, ye, tab_le, le, &leosTemp);
    if(status != ROOT_SUCCESS) {
      fprintf(stderr,
        "[EOS_SC_fill]: Failed to root find table!\n"
        "\trho      = %e\n"
        "\tu        = %e\n"
        "\te        = %e\n"
        "\tye       = %g\n"
        "\tlrho     = %g\n"
        "\tle       = %g\n"
        "\tlTguess  = %g\n"
        "\tleosTemp = %g\n",
        rho / GV::RHO_unit, u / GV::U_unit, e, ye, lrho, le, lTguess, leosTemp);
      exit(1); // TODO: Handle this more gracefully
    }

  if constexpr (SC_THROTTLE_CS){
    double cs2 = EOS_SC_sound_speed(lrho, lT, ye);
    if(cs2 >= 1.0) {
      if constexpr (SC_DEBUG){
        fprintf(stderr,
          "[EOS_SC_fill]: Warning! Sound speed superluminal!\n"
          "\tcs2   = %e\n"
          "\tlT    = %e\n"
          "\tlrho  = %e\n"
          "\tye    = %e\n"
          "\tNow throttling.\n",
          (cs2 / (C_LIGHT_CGS * C_LIGHT_CGS)), leosTemp, lrho, ye);
      }
      leosTemp = tab_lT_min;
      le = EOS_SC_interp(lrho, leosTemp, ye, tab_le);
      e = le2e(le);
      u = rho * e;
      assert(false && "uIn not found");
      //GV::uIn[UU] = u / GV::U_unit;
    }
  }

    eos[EOS_LRHO] = lrho;
    eos[EOS_LT] = leosTemp;
    eos[EOS_YE] = ye;

    if(isnan(leosTemp) == true)
      std::cout << "Particle id after leosTemp :      " << b.id() << std::endl;

    return;
  }

  static double
  EOS_SC_pressure_rho0_u(double lrho, double lT, double ye) {
    const double lP = EOS_SC_interp(lrho, lT, ye, tab_lP);
    return pow(10., lP) / GV::U_unit;
  }

  static double
  EOS_SC_specific_enthalpy_rho0_u(double lrho, double lT, double ye) {
    const double hm1 = EOS_SC_interp(lrho, lT, ye, tab_hm1);
    const double h_cgs = hm1 + C_LIGHT_CGS * C_LIGHT_CGS;
    const double h = h_cgs / (C_LIGHT_CGS * C_LIGHT_CGS);
    return h;
  }

  static double
  EOS_SC_sound_speed(double lrho, double lT, double ye) {
    double cs2 = EOS_SC_interp(lrho, lT, ye, tab_cs2);
    return sqrt(cs2);
  }

  static double
  EOS_SC_temperature(double lT) {
    // temperature is in MeV to start, which is a fine code unit
    // convert MeV to K
    return pow(10., lT) * MEV / KBOL; // / GV::TEMP_unit;
  }

  static double
  EOS_SC_entropy(double lrho, double lT, double ye) {
    double ent = EOS_SC_interp(lrho, lT, ye, tab_ent);
    return ent;
  }

  static double
  EOS_SC_gamma(double lrho, double lT, double ye) {
    const double cs2 = EOS_SC_interp(lrho, lT, ye, tab_cs2);
    const double lP = EOS_SC_interp(lrho, lT, ye, tab_lP);
    const double lgamma = lrho - lP;
    double gamma = cs2 * pow(10., lgamma);
    if(gamma < 1.1)
      gamma = 1.1;
    if(gamma > 2.0)
      gamma = 2.0;
    return gamma;
  }

  static double
  EOS_SC_get_u_of_T(double rho, double T, double ye) {
    if(T < tab_T_min)
      T = tab_T_min;
    if(T > tab_T_max)
      T = tab_T_max;
    const double lrho = catch_rho(rho);
    const double lT = catch_lT(log10(T));
    ye = catch_ye(ye);
    const double le = EOS_SC_interp(lrho, lT, ye, tab_le);
    const double e = le2e(le);
    const double u = e;
    return u / GV::U_unit;
  }

  static double
  EOS_SC_pressure_rho0_w(double rho, double w, double ye, double * lTold) {
    double lTguess = *lTold;
    double leosTemp;
    lTguess = catch_lT(lTguess);

    // Subtract off rest energy
    w -= rho;

    // convert to CGS
    rho *= GV::RHO_unit;
    w *= GV::U_unit;

    // dont' fall off the table
    double lrho = catch_rho(rho);
    double lw = catch_w(w);
    ye = catch_ye(ye);

    // force w back onto the table for chosen ilrho and iY
    lw = catch_var_2d(lrho, ye, lw, tab_lwmrho_min_2d, tab_lwmrho_max_2d);
  if constexpr (SC_MONOTONE_SAFE){
    lw = catch_var_2d_monotone(lrho, ye, lw, tab_lwmrho);
  }

  #if 0 // HL : again turn off
      // crash and die if something went wrong here
  if constexpr (SC_DEBUG){
      if (isnan(lw)) {
        fprintf(stderr,"[EOS_SC_Pressure_rho0_w %d]: NAN detected!\n",
          mpi_io_proc());
        fprintf(stderr,"rho     = %.10e\n",rho);
        fprintf(stderr,"w       = %.10e\n",w);
        fprintf(stderr,"lrho    = %.10f\n",lrho);
        fprintf(stderr,"lw      = %.10f\n",lw);
        fprintf(stderr,"ye      = %.10f\n",ye);
        fprintf(stderr,"lTguess = %.10f\n",lTguess);
        exit(1);
      }
  }
  #endif

    int status = find_lT(lrho, lTguess, ye, tab_lwmrho, lw, &leosTemp);

    if(status != ROOT_SUCCESS) {
      /*
      fprintf(stderr,
        "[EOS_SC_pressure_rho0_w %d]: "
        "Failed to root find table\n"
        "\trho      = %g\n"
        "\tlrho     = %g\n"
        "\tlTguess  = %g\n"
        "\tye       = %g\n"
        "\tleosTemp = %g\n"
        "\tlwmrho   = %g\n",
        mpi_io_proc(),rho,
        lrho, lTguess, ye, leosTemp,
        lw);
          */
      fprintf(stderr,
        "Failed to root find table\n"
        "\trho      = %g\n"
        "\tlrho     = %g\n"
        "\tlTguess  = %g\n"
        "\tye       = %g\n"
        "\tleosTemp = %g\n"
        "\tlwmrho   = %g\n",
        rho, lrho, lTguess, ye, leosTemp, lw);

      fprintf(stderr, "tab_lwwmrho = \n");
      temp_map(lrho, ye, tab_lwmrho);
      exit(1); // TODO: Handle this more gracefully
    }

    const double log_press = EOS_SC_interp(lrho, leosTemp, ye, tab_lP);
    double press = pow(10., log_press);

    // back into code units
    press /= GV::U_unit;

  if constexpr (SC_DEBUG){
    if(isnan(press)) {
      // TODO: handle this more gracefully.
      fprintf(stderr, "press from enthalpy = NaN.\n");
      exit(1);
    }
  }

    *lTold = leosTemp;
    return press;
  }

  static double
  EOS_SC_u_press(double press, double rho, double ye, double * lTold) {
    double lTguess = *lTold;
    double leosTemp; //, lrho, lp;
    lTguess = catch_lT(lTguess);

    // Units to CGS
    press *= GV::U_unit;
    rho *= GV::RHO_unit;

    const double lrho = catch_rho(rho);
    double lp = catch_press(press);
    ye = catch_ye(ye);
    lp = catch_var_2d(lrho, ye, lp, tab_lP_min_2d, tab_lP_max_2d);

  #if 0
      // crash and die if something went wrong here
  if constexpr (SC_DEBUG){
      if (isnan(lp)) {
        fprintf(stderr,"[EOS_SC_Pressure_rho0_w %d]: NAN detected!\n",
          mpi_io_proc());
        fprintf(stderr,"rho     = %.10e\n",rho);
        fprintf(stderr,"press   = %.10e\n",press);
        fprintf(stderr,"lrho    = %.10f\n",lrho);
        fprintf(stderr,"lP      = %.10f\n",lp);
        fprintf(stderr,"ye      = %.10f\n",ye);
        fprintf(stderr,"lTguess = %.10f\n",lTguess);
        exit(1);
      }
  }
  #endif

    int status = find_lT(lrho, lTguess, ye, tab_lP, lp, &leosTemp);

    if(status != ROOT_SUCCESS) {
      /*
        fprintf(stderr,"[EOS_SC_u_press]: "
          "Failed to root find table %g %g %g %g   %g %d\n",
                lrho, lTguess, ye, leosTemp, rho, mpi_io_proc());
      */
      fprintf(stderr,
        "[EOS_SC_u_press]: "
        "Failed to root find table %g %g %g %g %g\n",
        lrho, lTguess, ye, leosTemp, rho);
      exit(1); // TODO: Handle this more gracefully
    }
    const double le = EOS_SC_interp(lrho, leosTemp, ye, tab_le);
    const double e = le2e(le);
    double u = e;

    // to code units
    u /= GV::U_unit;

  if constexpr (SC_DEBUG){
    if(isnan(u)) {
      // TODO: handle this more gracefully.
      fprintf(stderr, "u from press = NaN.\n");
      fprintf(stderr, "press = %f\n", press);
      fprintf(stderr, "rho = %f\n", rho);
      fprintf(stderr, "ye = %f\n", ye);
      exit(1);
    }
  }

    *lTold = leosTemp;
    return u;
  }

  static void
  EOS_SC_mass_fractions(double Xi[NUM_MASS_FRACTIONS], const double * extra) {
    double lrho = extra[EOS_LRHO];
    double lT = extra[EOS_LT];
    double Ye = extra[EOS_YE];
    Xi[MF_XA] = EOS_SC_interp(lrho, lT, Ye, tab_Xa);
    Xi[MF_XH] = EOS_SC_interp(lrho, lT, Ye, tab_Xh);
    Xi[MF_XN] = EOS_SC_interp(lrho, lT, Ye, tab_Xn);
    Xi[MF_XP] = EOS_SC_interp(lrho, lT, Ye, tab_Xp);
  }

  static void
  EOS_SC_avg_ions(double * Abar, double * Zbar, const double * extra) {
    double lrho = extra[EOS_LRHO];
    double lT = extra[EOS_LT];
    double Ye = extra[EOS_YE];
    *Abar = EOS_SC_interp(lrho, lT, Ye, tab_Abar);
    *Zbar = EOS_SC_interp(lrho, lT, Ye, tab_Zbar);
  }

  static void
  EOS_SC_set_floors(double scale,
    double rho,
    double u,
    double ye,
    double bsq,
    double * rhoflr,
    double * uflr) {
    *rhoflr = EOS_SC_rho_floor(scale, bsq);
    *uflr = EOS_SC_u_floor(scale, bsq, ye);
    *rhoflr = std::max(*rhoflr, u / UORHOMAX);
  }

  static double
  EOS_SC_rho_floor(double scale, double bsq) {
    double rhoflr = RHOMIN * scale;
    double rhominlimit = RHOMINLIMIT;
    rhoflr = std::max(rhoflr, rhominlimit);
    rhoflr = std::max(rhoflr, bsq / BSQORHOMAX);
    return rhoflr;
  }

  static double
  EOS_SC_get_min_lrho() {
    return tab_lrho_min;
  }

  static double
  EOS_SC_get_min_rho() {
    double delrho = tab_lrho_max - tab_lrho_min;
    double lrho_min = tab_lrho_min + 0.01 * delrho;
    double rho_min_cgs = pow(10., lrho_min);
    double rho_min = rho_min_cgs / GV::RHO_unit;
    return rho_min;
  }

  static double
  EOS_SC_get_min_lT() {
    return tab_lT_min;
  }

  static double
  EOS_SC_get_minu(double rho, double ye) {
    double lrho = catch_rho(rho);
    double lT = tab_lT_min;
    ye = catch_ye(ye);

    double le = EOS_SC_interp(lrho, lT, ye, tab_le);
    double e = le2e(le);
    double u = e;
    return u / GV::U_unit;
  }

  static double
  EOS_SC_u_floor(double scale, double bsq, double ye) {
    // return -INFINITY;
    double rhoflr = EOS_SC_rho_floor(scale, bsq) * GV::RHO_unit;
    double minu = EOS_SC_get_minu(rhoflr, ye);
    if((bsq / BSQOUMAX) > fabs(minu)) { // Good idea?
      minu = bsq / BSQOUMAX;
    }
    return minu;
  }

  static void
  EOS_SC_get_polytrope(double lrho,
    double lT,
    double ye,
    double * poly_K,
    double * poly_gamma) {
    lrho = catch_lrho(lrho);
    lT = catch_lT(lT);
    ye = catch_ye(ye);
    double K = EOS_SC_interp(lrho, lT, ye, tab_poly_K);
    double Gam = EOS_SC_interp(lrho, lT, ye, tab_poly_gamma);
    double K_unit = GV::U_unit / pow(GV::RHO_unit, Gam);
    *poly_K = K / K_unit;
    *poly_gamma = Gam;
  }
  //----------------------------------------------------------------------

  // Root-finding
  // ----------------------------------------------------------------------
  static int
  find_lT(const double lrho,
    double lTguess,
    const double ye,
    double * tab,
    const double val,
    double * lT) {
    struct of_lT_params p;
    p.lrho = lrho;
    p.ye = ye;
    p.tab = tab;
    int status = find_root(&lT_f, &p, val, lTguess, tab_lT_min,
      tab_lT_max - 2. * TABLE_TOL / 100., TABLE_TOL, TABLE_FTOL, lT);
    return status;
  }

  static int
  find_adiabat_0d(double lrho, double lTguess, double ye, double s, double * lT) {
    lrho = catch_lrho(lrho);
    lTguess = catch_lT(lTguess);
    ye = catch_ye(ye);
    s = catch_s(s);

    struct of_lT_params p;
    p.lrho = lrho;
    p.ye = ye;
    p.tab = tab_ent;
    int status = find_root(&lT_f, &p, s, lTguess, tab_lT_min,
      tab_lT_max - 2. * TABLE_TOL / 100., TABLE_TOL, TABLE_FTOL, lT);
    return status;
  }

  static double
  EOS_SC_hm1_min_adiabat(const struct of_adiabat * a) {
    return a->hm1_min / (C_LIGHT_CGS * C_LIGHT_CGS);
  }

  static int
  EOS_SC_find_adiabat_1d(double s,
    double ye,
    double lrho_min,
    double lrho_max,
    struct of_adiabat * a) {
    s = catch_s(s);
    ye = catch_ye(ye);
    lrho_min = catch_lrho(lrho_min);
    lrho_max = catch_lrho(lrho_max);

    int ilrho_min = find_index(lrho_min, tab_lrho, Nrho);
    int ilrho_max = find_index(lrho_max, tab_lrho, Nrho);
    double slope = tab_dlT / tab_dlrho;
    double * lT_of_rho = safe_malloc<double>(Nrho);

    double hm1_min = INFINITY;
    double hm1_max = -INFINITY;

    for(int i = ilrho_min; i < ilrho_max; i++) {
      double lTguess = slope * tab_lrho[i];
      double lrho = tab_lrho[i];
      int status = find_adiabat_0d(lrho, lTguess, ye, s, &(lT_of_rho[i]));
      double hm1 = EOS_SC_interp(lrho, lT_of_rho[i], ye, tab_hm1);
      if(hm1 < hm1_min)
        hm1_min = hm1;
      if(hm1 > hm1_max)
        hm1_max = hm1;
      if(status != ROOT_SUCCESS)
        return ROOT_FAIL;
    }

    a->s = s;
    a->ye = ye;
    a->lrho_min = lrho_min;
    a->lrho_max = lrho_max;
    a->imin = ilrho_min;
    a->imax = ilrho_max;
    a->hm1_min = hm1_min;
    a->hm1_max = hm1_max;
    a->lT = lT_of_rho;

    return ROOT_SUCCESS;
  }

  static void
  EOS_SC_print_adiabat(const struct of_adiabat * a) {
    double ye = a->ye;
    fprintf(stdout, "\n");
    fprintf(stdout, "Adiabat =\n");
    fprintf(stdout, "#----------\n");
    fprintf(stdout, "#lrho\tlT\ts\n");
    for(int i = a->imin; i < a->imax; i++) {
      double lrho = tab_lrho[i];
      double lT = a->lT[i];
      fprintf(
        stdout, "%e\t%e\t%e\n", lrho, lT, EOS_SC_interp(lrho, lT, ye, tab_ent));
    }
    fprintf(stdout, "#----------\n");
    fprintf(stdout, "\n");
  }

  static void
  EOS_SC_adiabat_free(struct of_adiabat * a) {
    free(a->lT);
  }

  static void
  EOS_SC_isoentropy_hm1(double hm1,
    const struct of_adiabat * a,
    double * lrho_guess,
    double * rho,
    double * u) {
    hm1 = catch_hm1(hm1 * C_LIGHT_CGS * C_LIGHT_CGS);
    *lrho_guess = catch_lrho(*lrho_guess);
    double s = catch_s(a->s);
    double ye = catch_ye(a->ye);
    double lrho;

    // catch hm1 more carefully
    double hm1_min = a->hm1_min;
    if(hm1 < hm1_min)
      hm1 = hm1_min + TABLE_TOL / 100.;

    // root finding
    struct of_lT_adiabat_params p = {tab_hm1, a};
    int status = find_root(&lT_f_adiabat, &p, hm1, *lrho_guess, a->lrho_min,
      a->lrho_max, TABLE_TOL, TABLE_FTOL, &lrho);
    if(status != ROOT_SUCCESS) {
      fprintf(stderr,
        "[EOS_SC_isoentropy_hm1]: Failed to find root!\n"
        "\thm1        = %e\n"
        "\thm1_min    = %e\n"
        "\ts          = %e\n"
        "\tye         = %e\n"
        "\tlrho_guess = %e\n",
        hm1, hm1_min, s, ye, *lrho_guess);

      printf("ADIABAT MAP\n");
      printf("rho\thm1\n");
      for(int i = a->imin; i < a->imax; i++) {
        printf("%e\t%e\n", tab_lrho[i], lT_f_adiabat(tab_lrho[i], &p));
      }
      exit(1);
    }

    // and get what we care about
    double lT = interp_1d(
      lrho, a->lrho_min, a->lrho_max, a->imin, a->imax, tab_lrho, a->lT);
    double le = EOS_SC_interp(lrho, lT, ye, tab_le);
    double e = le2e(le);
    double rho_cgs = pow(10., lrho);
    double u_cgs = rho_cgs * e;
    // and to cgs
    *rho = rho_cgs / GV::RHO_unit;
    *u = u_cgs / GV::U_unit;
    // save initial guesses
    *lrho_guess = lrho;

    return;
  }
  // ----------------------------------------------------------------------

  // Interpolation
  // ----------------------------------------------------------------------
  static double
  EOS_SC_interp(const double lrho,
    const double lT,
    const double Ye,
    double * tab) {
    // We don't want our indices to fall off the table,
    // but we want to extrapolate appropriately
    double lT2 = catch_lT(lT);
    double Ye2 = catch_ye(Ye);
    double lrho2 = catch_lrho(lrho);

    // indices
    const int irho = (lrho2 - tab_lrho_min) / tab_dlrho;
    const int iT = (lT2 - tab_lT_min) / tab_dlT;
    const int iY = (Ye2 - tab_Ye_min) / tab_dYe;

    // assumes evenly spaced table
    const double delrho = (lrho - tab_lrho[irho]) / tab_dlrho;
    const double delT = (lT - tab_lT[iT]) / tab_dlT;
    const double delY = (Ye - tab_Ye[iY]) / tab_dYe;

    // one-sided, trilinear interpolation
    return ((1 - delrho) * (1 - delT) * (1 - delY) * tab[EOS_ELEM(irho, iT, iY)] +
            delrho * (1 - delT) * (1 - delY) * tab[EOS_ELEM(irho + 1, iT, iY)] +
            (1 - delrho) * delT * (1 - delY) * tab[EOS_ELEM(irho, iT + 1, iY)] +
            delrho * delT * (1 - delY) * tab[EOS_ELEM(irho + 1, iT + 1, iY)] +
            (1 - delrho) * (1 - delT) * delY * tab[EOS_ELEM(irho, iT, iY + 1)] +
            delrho * (1 - delT) * delY * tab[EOS_ELEM(irho + 1, iT, iY + 1)] +
            (1 - delrho) * delT * delY * tab[EOS_ELEM(irho, iT + 1, iY + 1)] +
            delrho * delT * delY * tab[EOS_ELEM(irho + 1, iT + 1, iY + 1)]);
  }

  static double
  interp_2d(const double lrho, const double Ye, const double * tab_2d) {
    // indices
    const int irho = (lrho - tab_lrho_min) / tab_dlrho;
    const int iY = (Ye - tab_Ye_min) / tab_dYe;

    // assumes evenly spaced table
    const double delrho = (lrho - tab_lrho[irho]) / tab_dlrho;
    const double delY = (Ye - tab_Ye[iY]) / tab_dYe;

    // one-sided, bilinear interpolation
    return ((1 - delrho) * (1 - delY) * tab_2d[MMA_ELEM(irho, iY)] +
            delrho * (1 - delY) * tab_2d[MMA_ELEM(irho + 1, iY)] +
            (1 - delrho) * delY * tab_2d[MMA_ELEM(irho, iY + 1)] +
            delrho * delY * tab_2d[MMA_ELEM(irho + 1, iY + 1)]);
  }
  // ----------------------------------------------------------------------

  // Error functions
  // ----------------------------------------------------------------------
  static double
  lT_f(const double lT, const void * params) {
    struct of_lT_params * p = (struct of_lT_params *)params;
    double * tab = p->tab;
    const double lrho = p->lrho;
    const double ye = p->ye;
    return EOS_SC_interp(lrho, lT, ye, tab);
  }

  static double
  lT_f_adiabat(const double lrho, const void * params) {
    struct of_lT_adiabat_params * p = (struct of_lT_adiabat_params *)params;
    const struct of_adiabat * a = p->a;
    double * tab = p->tab;

    const double ye = a->ye;
    const double lT = interp_1d(
      lrho, a->lrho_min, a->lrho_max, a->imin, a->imax, tab_lrho, a->lT);
  if constexpr (SC_DEBUG){
    if(isnan(lrho) || isnan(lT)) {
      fprintf(stderr,
        "[lT_f_adiabat]: NaN detected!\n"
        "\tlrho      = %f\n"
        "\tlT        = %f\n"
        "\tlrho_min  = %f\n"
        "\tlrho_max = %f\n"
        "\timin      = %d\n"
        "\timax      = %d\n"
        "\n",
        lrho, lT, a->lrho_min, a->lrho_max, a->imin, a->imax);
      exit(1);
    }
  }
    return EOS_SC_interp(lrho, lT, ye, tab);
  }

  static void
  EOS_SC_get_bounds(struct of_tablebounds * b) {
    b->Nrho = Nrho;
    b->NT = NT;
    b->NYe = NYe;

    b->lrho_min = tab_lrho_min;
    b->lrho_max = tab_lrho_max;
    b->dlrho = tab_dlrho;

    b->lT_min = tab_lT_min;
    b->lT_max = tab_lT_max;
    b->dlT = tab_dlT;

    b->Ye_min = tab_Ye_min;
    b->Ye_max = tab_Ye_max;
    b->dYe = tab_dYe;
  }
  // ----------------------------------------------------------------------

  // Utilities
  // ----------------------------------------------------------------------
  static void
  fill_min_1d(double * tab_min_1d, double * tab) {
    int elem;
    for(int iY = 0; iY < NYe; iY++) {
      double min = INFINITY;
      for(int iT = 0; iT < NT; iT++) {
        for(int irho = 0; irho < Nrho; irho++) {
          elem = EOS_ELEM(irho, iT, iY);
          if(tab[elem] < min)
            min = tab[elem];
        }
      }
      tab_min_1d[iY] = min;
    }
  }

  // Reduces the 3d tables to 2d max/min tables
  // where the maximum is taken over lT for each pair of lrho,ye values
  static void
  fill_max_min_2d(double * tab_max_2d, double * tab_min_2d, double * tab) {
    int elem, elem2d;
    for(int irho = 0; irho < Nrho; irho++) {
      for(int iY = 0; iY < NYe; iY++) {
        elem2d = MMA_ELEM(irho, iY);
        double max = -INFINITY;
        double min = INFINITY;
        for(int iT = 0; iT < NT; iT++) {
          elem = EOS_ELEM(irho, iT, iY);
          if(tab[elem] > max)
            max = tab[elem];
          if(tab[elem] < min)
            min = tab[elem];
        }
        tab_max_2d[elem2d] = max;
        tab_min_2d[elem2d] = min;
      }
    }
  }

  static double
  catch_var_2d(const double lrho,
    const double Ye,
    const double var,
    const double * tab_min_2d,
    const double * tab_max_2d) {
  if constexpr (SC_DEBUG){
    {
      if(isnan(lrho) || isnan(Ye) || isnan(var)) {
        fprintf(stderr,
          "[EOS_SC::catch_var_2d]: NaN detected!\n"
          "\tlrho = %e\n"
          "\tYe   = %e\n"
          "\tvar  = %e\n",
          lrho, Ye, var);
      }
    }
  }

    const double varmin = interp_2d(lrho, Ye, tab_min_2d);
    if(var <= varmin)
      return varmin + TABLE_TOL / 100.;
    const double varmax = interp_2d(lrho, Ye, tab_max_2d);
    if(var >= varmax)
      return varmax - TABLE_TOL / 100.;
    return var;
  }

  static double
  catch_var_2d_monotone(const double lrho,
    const double Ye,
    const double var,
    double * tab) {
    const double varmin = EOS_SC_interp(lrho, tab_lT_min, Ye, tab);
    if(var <= varmin)
      return varmin + TABLE_TOL / 100.;
    const double varmax = EOS_SC_interp(lrho, tab_lT_max, Ye, tab);
    if(var >= varmax)
      return varmax - TABLE_TOL / 100.;
    return var;
  }

  static void
  temp_map(double lrho, double Ye, const double * tab) {
    // Don't fall off table
    if(lrho < tab_lrho_min)
      lrho = tab_lrho_min;
    if(lrho >= tab_lrho_max)
      lrho = tab_lrho_max - TABLE_TOL / 100.;
    if(Ye < tab_Ye_min)
      Ye = tab_Ye_min;
    if(Ye >= tab_Ye_max)
      Ye = tab_Ye_max - TABLE_TOL / 100.;
    // Indices
    const int irho = (lrho - tab_lrho_min) / tab_dlrho;
    const int iY = (Ye - tab_Ye_min) / tab_dYe;
    for(int iT = 0; iT < NT; iT++) {
      fprintf(stderr, "%d\t%.10f\t%.10f\n", iT, tab_lT[iT],
        tab[EOS_ELEM(irho, iT, iY)]);
    }
  }

  static double
  le2e(const double le) {
    return pow(10., le) - energy_shift;
  }

  static double
  e2le(const double e) {
    return log10(e + energy_shift);
  }

  /*
  static double lw2w(const double lw)
  {
    return pow(10., lw) - enthalpy_shift;
  }
  */

  static double
  w2lw(const double w) {
    return log10(w + enthalpy_shift);
  }

  static double
  catch_rho(const double rho) {
    // dont' fall off the table
    if(rho < tab_rho_min) {
      return tab_lrho_min;
    }
    else if(rho >= tab_rho_max) {
      return (tab_lrho_max - TABLE_TOL / 100.);
    }
    else {
      return log10(rho);
    }
  }

  static double
  catch_lrho(const double lrho) {
    if(lrho < tab_lrho_min) {
      return tab_lrho_min;
    }
    else if(lrho > tab_lrho_max) {
      return tab_lrho_max - TABLE_TOL / 100.;
    }
    else {
      return lrho;
    }
  }

  static double
  catch_e(const double e) {
    if(e < tab_e_min) {
      return tab_le_min;
    }
    else if(e >= tab_e_max) {
      return tab_le_max - TABLE_TOL / 100.;
    }
    else {
      return e2le(e);
    }
  }

  static double
  catch_press(const double press) {
    if(press < tab_P_min) {
      return tab_lP_min;
    }
    else if(press >= tab_P_max) {
      return tab_lP_max - TABLE_TOL / 100.;
    }
    else {
      return log10(press);
    }
  }

  static double
  catch_w(const double w) {
    if(w < tab_wmrho_min) {
      return tab_lwmrho_min;
    }
    else if(w >= tab_wmrho_max) {
      return tab_lwmrho_max - TABLE_TOL / 100.;
    }
    else {
      return w2lw(w);
    }
  }

  /*
  static double catch_h(const double h)
  {
    if ( h < tab_hm1_min ) {
      return tab_hm1_min;
    } else if ( h >= tab_hm1_max ) {
      return tab_hm1_max - TABLE_TOL/100.;
    } else {
      return h;
    }
  }
  */

  static double
  catch_ye(const double ye) {
    if(ye < tab_Ye_min) {
      return tab_Ye_min;
    }
    else if(ye >= tab_Ye_max) {
      return tab_Ye_max - TABLE_TOL / 100.;
    }
    else {
      return ye;
    }
  }

  /*
  static double catch_temp(const double temp)
  {
    if ( temp < tab_T_min ) {
      return tab_lT_min;
    } else if ( temp >= tab_T_max ) {
      return tab_lT_max - TABLE_TOL/100.;
    } else {
      return log10(temp);
    }
  }
  */

  static double
  catch_lT(const double lT) {
    if(lT < tab_lT_min) {
      return tab_lT_min;
    }
    else if(lT >= tab_lT_max) {
      return tab_lT_max - TABLE_TOL / 100.;
      // } else if (isnan(lT)) {
      //   return 0.5*(tab_lT_min + tab_lT_max);
    }
    else {
      return lT;
    }
  }

  static double
  catch_s(const double s) {
    if(s < tab_ent_min) {
      return tab_ent_min;
    }
    else if(s >= tab_ent_max) {
      return tab_ent_max - TABLE_TOL / 100.;
    }
    else {
      return s;
    }
  }

  static double
  catch_hm1(const double hm1) {
    if(hm1 < tab_hm1_min) {
      return tab_hm1_min;
    }
    else if(hm1 >= tab_hm1_max) {
      return tab_hm1_max - TABLE_TOL / 100.;
    }
    else {
      return hm1;
    }
  }
  // ----------------------------------------------------------------------
  //
  //
  // from eos.c. This contains all main functionalities after above SC
  // readers

  static double
  EOS_bad_eos_error() {
    fprintf(stderr, "ERROR! UNKNOWN EOS TYPE!\n");
    exit(1);
  }

  /*******************************************************************************
        Wrappers
  *******************************************************************************/

  // Getting pressure from rho and u
  // eos_cache saves unevolved additional variables.
  static double
  EOS_pressure_rho0_u(body & b) {
    // Call EOS_SC_fill that put neccessary cache data.
    // Log(T) is calculated and save from this function
    double eos_cache[3];
    EOS_SC_fill(b, eos_cache);

    // Particles data
    double rho = b.getDensity();
    double u = b.getInternalenergy();
    double ye = b.getElectronfraction();
    double lrho = eos_cache[EOS_LRHO];
    double lT = eos_cache[EOS_LT];

    double press;
    double rho_poly_thresh = EOS_SC_get_min_rho();
    if(rho < rho_poly_thresh) {
      // double K,Gam;
      // lrho = EOS_SC_get_min_lrho();
      // lT = EOS_SC_get_min_lT();
      // EOS_SC_get_polytrope(lrho, lT, ye, &K, &Gam);
      // press = EOS_Poly_pressure_rho0_u(rho,u,K,Gam);
      // // double press_min =
      // EOS_SC_pressure_rho0_u(log(rho_poly_thresh),log(b.getTemperature()*KBOL/MEV),ye);//*rho/rho_poly;
      // // double r_min =
      // density_profiles::r_from_rho_grid_input_file(rho_poly_thresh);
      // // double press_min = density_profiles::p_from_input_file
      double lrho_thresh = EOS_SC_get_min_lrho();
      double press_min = EOS_SC_pressure_rho0_u(lrho_thresh, lT, ye);
      press = press_min * pow(rho/rho_poly_thresh, param::gamma_poly_thresh);
    }
    else {
      press = EOS_SC_pressure_rho0_u(lrho, lT, ye);
    }
  #if 0
        EOS_bad_eos_error();
  #endif
    return press;
  }

  static double
  EOS_enthalpy_rho0_u(double rho, double u, const double * extra) {
    double enth;
    double lrho = extra[EOS_LRHO];
    double lT = extra[EOS_LT];
    double ye = extra[EOS_YE];
    double rho_poly_thresh = EOS_SC_get_min_rho();
    if(rho < rho_poly_thresh) {
      double K, Gam;
      lrho = EOS_SC_get_min_lrho();
      lT = EOS_SC_get_min_lT();
      EOS_SC_get_polytrope(lrho, lT, ye, &K, &Gam);
      enth = EOS_Poly_enthalpy_rho0_u(rho, std::max(u, 0.0), K, Gam);
    }
    else {
      double h = EOS_SC_specific_enthalpy_rho0_u(lrho, lT, ye);
      enth = h * rho;
    }
  #if 0
        EOS_bad_eos_error();
  #endif
    return enth;
  }

  static double
  EOS_entropy_rho0_u(double rho, double u, const double * extra) {
    double ent;
    double lrho = extra[EOS_LRHO];
    double lT = extra[EOS_LT];
    double ye = extra[EOS_YE];
    double rho_poly_thresh = EOS_SC_get_min_rho();
    if(rho < rho_poly_thresh) {
      double K, Gam;
      lrho = EOS_SC_get_min_lrho();
      lT = EOS_SC_get_min_lT();
      EOS_SC_get_polytrope(lrho, lT, ye, &K, &Gam);
      ent = EOS_Poly_entropy_rho0_u(rho, u, K, Gam);
    }
    else {
      ent = EOS_SC_entropy(lrho, lT, ye);
    }
  #if 0
        EOS_bad_eos_error();
  #endif
    return ent;
  }

  // double rho, double u, const double* extra
  static double
  EOS_sound_speed_rho0_u(body & b) {

    // Call EOS_SC_fill to save eos_cache data
    double eos_cache[3];
    EOS_SC_fill(b, eos_cache);

    // particle data
    double rho = b.getDensity();
    double lrho = eos_cache[EOS_LRHO];
    double u = b.getInternalenergy();
    double ye = b.getElectronfraction();
    double lT = eos_cache[EOS_LT];

    double cs;
    double rho_poly_thresh = EOS_SC_get_min_rho();
    if(rho < rho_poly_thresh) {
      // double K,Gam;
      // lrho = EOS_SC_get_min_lrho();
      // lT = EOS_SC_get_min_lT();
      // EOS_SC_get_polytrope(lrho, lT, ye, &K, &Gam);
      // cs = EOS_Poly_sound_speed_rho0_u(rho,u,K,Gam);
      cs = sqrt(param::poly_gamma * b.getPressure() / b.getDensity());
    }
    else {
      cs = EOS_SC_sound_speed(lrho, lT, ye);
    }
  #if 0 // HL : will put correct conditional statement to call bad eos
        EOS_bad_eos_error();
  #endif
    return cs;
  }

  // double rho, double u, const double* extra
  static double
  EOS_temperature_sc(body & b) {

    // Call EOS_SC_fill to save eos_cache data
    double eos_cache[3];
    EOS_SC_fill(b, eos_cache);

    // particle data
    double rho = b.getDensity();
    double lrho = eos_cache[EOS_LRHO];
    double u = b.getInternalenergy();
    double ye = b.getElectronfraction();
    double lT = eos_cache[EOS_LT];

    // double cs;
    double temp;
    double rho_poly_thresh = EOS_SC_get_min_rho();
    if(rho < rho_poly_thresh) {
      //   double K,Gam;
      //   lrho = EOS_SC_get_min_lrho();
      //   lT = EOS_SC_get_min_lT();
      //   EOS_SC_get_polytrope(lrho, lT, ye, &K, &Gam);
      //   cs = EOS_Poly_sound_speed_rho0_u(rho,u,K,Gam);
      temp = b.getPressure() / (b.getDensity() * 8.314e7);
    }
    else {
      // cs = EOS_SC_sound_speed(lrho,lT,ye);
      temp = pow(10, lT) * MEV / KBOL;
    }
  #if 0 // HL : will put correct conditional statement to call bad eos
        EOS_bad_eos_error();
  #endif
    // return cs;
    return temp;
  }

  static void
  EOS_set_floors(double scale,
    double rho,
    double u,
    double bsq,
    double * rhoflr,
    double * uflr,
    const double * extra) {
    double ye = extra[EOS_YE];
    EOS_SC_set_floors(scale, rho, u, ye, bsq, rhoflr, uflr);
  #if 0
      EOS_bad_eos_error();
  #endif
  }

  static double
  EOS_adiabatic_constant(double rho, double u, const double * extra) {
    double cad;
    double gam = EOS_get_gamma(extra);
    cad = u * pow(rho, -gam);
  #if 0
      EOS_bad_eos_error();
  #endif
    return cad;
  }

  static double
  EOS_get_gamma(const double * extra) {
    double lrho = extra[EOS_LRHO];
    double lT = extra[EOS_LT];
    double ye = extra[EOS_YE];
    double gam = EOS_SC_gamma(lrho, lT, ye);
    return gam;
  #if 0
      EOS_bad_eos_error();
  #endif
  }

  static double
  EOS_temperature(double rho, double u, const double * extra) {
    double lT = extra[EOS_LT];
    return EOS_SC_temperature(lT);
  #if 0
      EOS_bad_eos_error();
  #endif
  }

  static double
  EOS_u_press(double press, double rho, double * extra) {
    double u;
    double ye = extra[EOS_YE];
    // double yedens = extra[EOS_YE];
    // double ye     = fabs(yedens) / (fabs(rho) + SMALL);
    double lTold = extra[EOS_LT];
    double rho_poly_thresh = EOS_SC_get_min_rho();
    if(rho < rho_poly_thresh) {
      u = EOS_SC_get_minu(rho, ye);
      lTold = EOS_SC_get_min_lT();
    }
    else {
      u = EOS_SC_u_press(press, rho, ye, &lTold);
    }
    extra[EOS_LT] = lTold;
  #if 0
      EOS_bad_eos_error();
  #endif // EOS
    return u;
  }

  static void
  EOS_SC_overwrite_cs2_with_table(const char * name) {

    std::ifstream infile(name);
    if(!infile.good()) {
      log_one(error) << "File " << name << " not found." << std::endl;
      MPI_Finalize();
      exit(-1);
    }

    hid_t file_id = H5Fopen(name, H5F_ACC_RDONLY, H5P_DEFAULT);

    hid_t dset_id = H5Dopen(file_id, "cs2", H5P_DEFAULT);
    int status =
      H5Dread(dset_id, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, tab_cs2);
    H5Dclose(dset_id);

    H5Fclose(file_id);

    int tab_size = NYe * NT * Nrho;
    tab_cs2_min = find_min(tab_cs2, tab_size);
    tab_cs2_max = find_max(tab_cs2, tab_size);
  }

  double
  EOS_temperature(body & b) {
    // Call EOS_SC_fill to save eos_cache data
    double eos_cache[3];
    EOS_SC_fill(b, eos_cache);

    double lT = eos_cache[EOS_LT];

    // Return temperature in CGS units
    return EOS_SC_temperature(lT);
  #if 0
      EOS_bad_eos_error();
  #endif
  }


  static void
  EOS_root_find_failure_test() {
    std::cout << "Reading EOS from file: " << param::eos_tab_file_path
              << std::endl;

    // Read in the table
    std::cout << "Entering function EOS_SC_init..." << std::endl;
    EOS_SC_init(param::eos_tab_file_path);

    // generate the randomly selected sample inputs
    const int num_samples = 50;
    std::vector<double> le_samples(std::pow(num_samples, 3));
    std::vector<double> lr_samples(num_samples);
    std::vector<double> Ye_samples(num_samples);
    std::vector<double> lT_samples(num_samples);
    std::vector<double> lT_roots(std::pow(num_samples, 3));
    double shift = 0.99; // how much of the table to sample?
    double lr_sample_width = (tab_lrho_max - tab_lrho_min);
    double Ye_sample_width = (tab_Ye_max - tab_Ye_min);
    double lT_sample_width = (tab_lT_max - tab_lT_min);
    for(int i = 0; i < num_samples; i++) {
      double prop = (double)i / (double)num_samples;
      // generate the le input
      // le_samples[i]=tab_le_min+(tab_le_max-tab_le_min)*prop;
      // le_samples[i]=catch_e(le2e(le_samples[i]));
      // generate the density input
      lr_samples[i] =
        tab_lrho_min + lr_sample_width * (0.5 * (1. - shift) + prop * shift);
      lr_samples[i] = catch_lrho(lr_samples[i]);
      // generate the Ye input
      Ye_samples[i] =
        tab_Ye_min + Ye_sample_width * (0.5 * (1. - shift) + prop * shift);
      Ye_samples[i] = catch_ye(Ye_samples[i]);
      // generate the lT input
      lT_samples[i] =
        tab_lT_min + lT_sample_width * (0.5 * (1. - shift) + prop * shift);
      lT_samples[i] = catch_lT(lT_samples[i]);
    }

    // test the root finder
    int idx = 0, // index for log(e) samples and log(T) roots
      num_fails = 0, // index to count incorrect roots
      num_degen = 0, // index to count root finder failures
      num_records = 0; // index to count number of records (fails+degen+bounds)
    double success_rel_err = 0.0, // relative error variable
      success_abs_err = 0.0; // absolute error variable
    double lT_guess =
      catch_lT(0.0); // the guess is zero in EOS_stellar_collapse.h
    std::vector<table_val>
      records; // vector to hold records (printed to txt file)
    table_val tv; // struct to hold relevant information

    // loop over all the samples
    for(int lT_idx = 0; lT_idx < num_samples; lT_idx++) {
      for(int lr_idx = 0; lr_idx < num_samples; lr_idx++) {
        for(int Ye_idx = 0; Ye_idx < num_samples; Ye_idx++) {

          // interpolate the energy from the log rho, log T, Ye sample point
          le_samples[idx] = EOS_SC_interp(
            lr_samples[lr_idx], lT_samples[lT_idx], Ye_samples[Ye_idx], tab_le);

          // root find lT from the log rho, Ye, and interpolated le values
          int status = find_lT(lr_samples[lr_idx], lT_guess, Ye_samples[Ye_idx],
            tab_le, le_samples[idx], &lT_roots[idx]);

          // calculate the relative error
          // for checking if correct root is found
          double abs_tmp = std::fabs(lT_samples[lT_idx] - lT_roots[idx]);
          double rel_tmp = abs_tmp / std::fabs(lT_samples[lT_idx]);
          if(abs_tmp > success_abs_err && status == ROOT_SUCCESS)
            success_abs_err = abs_tmp;
          if(rel_tmp > success_rel_err && status == ROOT_SUCCESS)
            success_rel_err = rel_tmp;

          // Found the wrong root
          // The root finder reported success
          // but there is a large difference between the original temp
          // and the root find temperature
          if(rel_tmp > 1e-6 && status == ROOT_SUCCESS) {
            num_degen++;
            num_records++;
            tv.lrho = lr_samples[lr_idx];
            tv.lT = lT_samples[lT_idx];
            tv.Ye = Ye_samples[Ye_idx];
            tv.le = le_samples[idx];
            tv.lT_root = lT_roots[idx];
            tv.rt = degen;
            tv.lr_idx = lr_idx;
            tv.Ye_idx = Ye_idx;
            records.push_back(tv);
          }

          // Could not find a root
          // secant method and bisection method failed
          if(status != ROOT_SUCCESS) {
            num_fails++;
            num_records++;
            tv.lrho = lr_samples[lr_idx];
            tv.lT = lT_samples[lT_idx];
            tv.Ye = Ye_samples[Ye_idx];
            tv.le = le_samples[idx];
            tv.lT_root = lT_roots[idx];
            tv.rt = fail;
            tv.lr_idx = lr_idx;
            tv.Ye_idx = Ye_idx;
            records.push_back(tv);
          }

          // record the boundary values for plotting purposes
          if(lT_idx == 0 || lT_idx == num_samples - 1 || lr_idx == 0 ||
            lr_idx == num_samples - 1 || Ye_idx == 0 ||
            Ye_idx == num_samples - 1) {
            num_records++;
            tv.lrho = lr_samples[lr_idx];
            tv.lT = lT_samples[lT_idx];
            tv.Ye = Ye_samples[Ye_idx];
            tv.le = le_samples[idx];
            tv.lT_root = lT_roots[idx];
            tv.rt = success_boundary;
            tv.lr_idx = lr_idx;
            tv.Ye_idx = Ye_idx;
            records.push_back(tv);
          }
          idx++; // icrement le_samples, lT_roots
        } // lrho_idx loop
      } // Ye_idx loop
    } // le_idx loop

    // open the output file (goes to build/test/test)
    std::ofstream outfile("root_find_failures.txt");
    if(outfile.is_open()) {
      std::cout << "Opened successfully" << std::endl;
    }
    else {
      std::cout << "Failed to open output file" << std::endl;
    }
    // write the records
    for(int i = 0; i < num_records; i++) {
      outfile << i << "," << records[i].lrho << "," << records[i].lT << ","
              << records[i].Ye << "," << records[i].le << ","
              << records[i].lT_root << "," << records[i].rt << ","
              << records[i].lr_idx << "," << records[i].Ye_idx << ","
              << num_samples << "," << num_records;
      // for a given lrho, Ye combination
      // data to plot le = func(lT)
      for(int j = 0; j < num_samples; j++) {
        double out_le =
          EOS_SC_interp(records[i].lrho, lT_samples[j], records[i].Ye, tab_le);
        outfile << "," << lT_samples[j] << "," << out_le;
      }
      outfile << std::endl;
    }
    outfile.close();
  }
};

#if eos_type == eos_stellar_collapse
  get_quantity_t eos_t<param::eos_stellar_collapse>::get_dpdrho_at_temp = nullptr;
#endif


// Init variable
int eos_t<param::eos_stellar_collapse>::Nrho = 0;
int eos_t<param::eos_stellar_collapse>::NT = 0;
int eos_t<param::eos_stellar_collapse>::NYe = 0;
int eos_t<param::eos_stellar_collapse>::Nrho_ye = 0;
int eos_t<param::eos_stellar_collapse>::NT_ye = 0;
int eos_t<param::eos_stellar_collapse>::NYe_ye = 0;

double * eos_t<param::eos_stellar_collapse>::tab_lrho = nullptr;
double * eos_t<param::eos_stellar_collapse>::tab_lT = nullptr;
double * eos_t<param::eos_stellar_collapse>::tab_Ye = nullptr;
double * eos_t<param::eos_stellar_collapse>::tab_lP = nullptr;
double * eos_t<param::eos_stellar_collapse>::tab_ent = nullptr;
double * eos_t<param::eos_stellar_collapse>::tab_dpderho = nullptr;
double * eos_t<param::eos_stellar_collapse>::tab_dpdrhoe = nullptr;
double * eos_t<param::eos_stellar_collapse>::tab_cs2 = nullptr;
double * eos_t<param::eos_stellar_collapse>::tab_le = nullptr;
double * eos_t<param::eos_stellar_collapse>::tab_Xa = nullptr;
double * eos_t<param::eos_stellar_collapse>::tab_Xh = nullptr;
double * eos_t<param::eos_stellar_collapse>::tab_Xn = nullptr;
double * eos_t<param::eos_stellar_collapse>::tab_Xp = nullptr;
double * eos_t<param::eos_stellar_collapse>::tab_Abar = nullptr;
double * eos_t<param::eos_stellar_collapse>::tab_Zbar = nullptr;
double * eos_t<param::eos_stellar_collapse>::tab_lwmrho = nullptr; // log enthalpy - rho, by volume
double * eos_t<param::eos_stellar_collapse>::tab_hm1 = nullptr; // enthalpy - 1, by mass
double * eos_t<param::eos_stellar_collapse>::tab_poly_gamma = nullptr; // Polytrope gamma
double * eos_t<param::eos_stellar_collapse>::tab_poly_K = nullptr; // polytrope K
double * eos_t<param::eos_stellar_collapse>::tab_rho = nullptr;
double * eos_t<param::eos_stellar_collapse>::tab_Yeye = nullptr;
double * eos_t<param::eos_stellar_collapse>::tab_T = nullptr;
double * eos_t<param::eos_stellar_collapse>::tab_dYedt = nullptr;
double * eos_t<param::eos_stellar_collapse>::tab_deweakdt = nullptr;

// min and max of wmrho given fixed ilrho and iY
double * eos_t<param::eos_stellar_collapse>::tab_le_min_2d = 0;
double * eos_t<param::eos_stellar_collapse>::tab_le_max_2d = 0;
double * eos_t<param::eos_stellar_collapse>::tab_lP_min_2d = 0;
double * eos_t<param::eos_stellar_collapse>::tab_lP_max_2d = 0;
double * eos_t<param::eos_stellar_collapse>::tab_lwmrho_min_2d = 0;
double * eos_t<param::eos_stellar_collapse>::tab_lwmrho_max_2d = 0;
double * eos_t<param::eos_stellar_collapse>::tab_hm1_min_1d = 0;

double eos_t<param::eos_stellar_collapse>::tab_lrho_min = 0;
double eos_t<param::eos_stellar_collapse>::tab_lrho_max = 0;
double eos_t<param::eos_stellar_collapse>::tab_rhoye_min = 0;
double eos_t<param::eos_stellar_collapse>::tab_rhoye_max = 0;
double eos_t<param::eos_stellar_collapse>::tab_lT_min = 0;
double eos_t<param::eos_stellar_collapse>::tab_lT_max = 0;
double eos_t<param::eos_stellar_collapse>::tab_Tye_min = 0;
double eos_t<param::eos_stellar_collapse>::tab_Tye_max = 0;
double eos_t<param::eos_stellar_collapse>::tab_Ye_min = 0;
double eos_t<param::eos_stellar_collapse>::tab_Ye_max = 0;
double eos_t<param::eos_stellar_collapse>::tab_Yeye_min = 0;
double eos_t<param::eos_stellar_collapse>::tab_Yeye_max = 0;
double eos_t<param::eos_stellar_collapse>::tab_dlrho = 0;
double eos_t<param::eos_stellar_collapse>::tab_dlT = 0;
double eos_t<param::eos_stellar_collapse>::tab_dYe = 0;
double eos_t<param::eos_stellar_collapse>::tab_drhoye = 0;
double eos_t<param::eos_stellar_collapse>::tab_dTye = 0;
double eos_t<param::eos_stellar_collapse>::tab_dYeye = 0;

double eos_t<param::eos_stellar_collapse>::tab_lP_min = 0;
double eos_t<param::eos_stellar_collapse>::tab_lP_max = 0;
double eos_t<param::eos_stellar_collapse>::tab_ent_min = 0;
double eos_t<param::eos_stellar_collapse>::tab_ent_max = 0;
double eos_t<param::eos_stellar_collapse>::tab_cs2_min = 0;
double eos_t<param::eos_stellar_collapse>::tab_cs2_max = 0;
double eos_t<param::eos_stellar_collapse>::tab_le_min = 0;
double eos_t<param::eos_stellar_collapse>::tab_le_max = 0;
double eos_t<param::eos_stellar_collapse>::tab_Xa_min = 0;
double eos_t<param::eos_stellar_collapse>::tab_Xa_max = 0;
double eos_t<param::eos_stellar_collapse>::tab_Xh_min = 0;
double eos_t<param::eos_stellar_collapse>::tab_Xh_max = 0;
double eos_t<param::eos_stellar_collapse>::tab_Xn_min = 0;
double eos_t<param::eos_stellar_collapse>::tab_Xn_max = 0;
double eos_t<param::eos_stellar_collapse>::tab_Xp_min = 0;
double eos_t<param::eos_stellar_collapse>::tab_Xp_max = 0;
double eos_t<param::eos_stellar_collapse>::tab_Abar_min = 0;
double eos_t<param::eos_stellar_collapse>::tab_Abar_max = 0;
double eos_t<param::eos_stellar_collapse>::tab_Zbar_min = 0;
double eos_t<param::eos_stellar_collapse>::tab_Zbar_max = 0;
double eos_t<param::eos_stellar_collapse>::tab_dpderho_min = 0;
double eos_t<param::eos_stellar_collapse>::tab_dpderho_max = 0;
double eos_t<param::eos_stellar_collapse>::tab_dpdrhoe_min = 0;
double eos_t<param::eos_stellar_collapse>::tab_dpdrhoe_max = 0;
double eos_t<param::eos_stellar_collapse>::tab_lwmrho_min = 0;
double eos_t<param::eos_stellar_collapse>::tab_lwmrho_max = 0;
double eos_t<param::eos_stellar_collapse>::tab_dYedt_min = 0;
double eos_t<param::eos_stellar_collapse>::tab_dYedt_max = 0;
double eos_t<param::eos_stellar_collapse>::tab_deweakdt_min = 0;
double eos_t<param::eos_stellar_collapse>::tab_deweakdt_max = 0;

double eos_t<param::eos_stellar_collapse>::tab_rho_min = 0;
double eos_t<param::eos_stellar_collapse>::tab_rho_max = 0;
double eos_t<param::eos_stellar_collapse>::tab_T_min = 0;
double eos_t<param::eos_stellar_collapse>::tab_T_max = 0;
double eos_t<param::eos_stellar_collapse>::tab_e_min = 0;
double eos_t<param::eos_stellar_collapse>::tab_e_max = 0;
double eos_t<param::eos_stellar_collapse>::tab_P_min = 0;
double eos_t<param::eos_stellar_collapse>::tab_P_max = 0;
double eos_t<param::eos_stellar_collapse>::tab_wmrho_min = 0;
double eos_t<param::eos_stellar_collapse>::tab_wmrho_max = 0;
double eos_t<param::eos_stellar_collapse>::tab_hm1_min = 0;
double eos_t<param::eos_stellar_collapse>::tab_hm1_max = 0;

double eos_t<param::eos_stellar_collapse>::pressure_min = 0;

double eos_t<param::eos_stellar_collapse>::energy_shift = 0;
double eos_t<param::eos_stellar_collapse>::enthalpy_shift = 0;

} // namespace eos
