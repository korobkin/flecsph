/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2018 Triad National Security, LLC
 * All rights reserved.
 * --------------------------------------------------------------------------~*/

/******************************************************************************
 *                                                                            *
 * EOS_HELM.h                                                                 *
 * TODO: CLEANUP, FIX MEMORY ALLOCATIONS                                      *
 *                                                                            *
 * IMPLEMENTS ROUTINES FOR READING EOS TABLES PROVIDED ON COCOCUBED.COM       *
 * ALSO INCLUDES FUNCTIONS FOR CALCULATING EOS VALUES                         *
 *                                                                            *
 ******************************************************************************/
#pragma once

#include "params.h"
#include "units.h"
#include <fstream>
#include "interp.h"
#define SQ(x) ((x) * (x))
#define CU(x) ((x) * (x) * (x))
#define QU(x) ((x) * (x) * (x) * (x))

#ifndef HELM_EOS_MAXITER
#define HELM_EOS_MAXITER 50
#endif // HELM_EOS_MAXITER
#ifndef HELM_EOS_EPS
#define HELM_EOS_EPS 1e-10
#endif // HELM_EOS_EPS

// Exclude Coulomb screening correction
#define HELM_EOS_NO_COULOMB_SCREENING
//#define IMAX 541
//#define JMAX 201
namespace eos {

namespace Detail {
    double constexpr
    expTaylor(double x, int n) {
        return n == 32
            ? x/32.
            : 1. + x/n*expTaylor(x, n + 1);
    }

    long constexpr
    exp10_int(int n) {
        return n == 0
            ? 1
            : 10*exp10_int(n - 1);
    }

}

/*
 * Constexpr version of the exponent-10
 */
double constexpr
exp10_constexpr(double x) {
    constexpr double ln10 = 2.30258509299404568401799145468436420760;
    return x > -1
        ?    Detail::expTaylor((x - (long)x)*ln10, 1)
           * Detail::exp10_int((long)(x))
        : 1./Detail::expTaylor(((long)x - x)*ln10, 1)
            /Detail::exp10_int((long)(-x));
}


using std::log; // to avoid name clash with flecsi::log

template<>
class eos_t<param::eos_helmholtz> {

  enum state_values {
    VALUE = 0,
    DRHO  = 1,
    DTEMP = 2,
    DABAR = 3,
    DZBAR = 4
  };

  // Hardcoded Helmholtz EOS table parameters
  // TODO: hardcoding bad, remove sometimes
  static constexpr int
    //tab_nrho = 271,
    //tab_ntemp = 101;
    tab_nrho = 541,
    tab_ntemp = 201,
    tab_extrapolation_margin_irho = 10,
    tab_extrapolation_margin_itemp = 10;

  static constexpr double
    tab_ltemp_min = 3.,
    tab_ltemp_max = 13.,
    tab_lrho_min  =-12.,
    tab_lrho_max  = 15.,
    tab_temp_min = exp10_constexpr(tab_ltemp_min),
    tab_temp_max = exp10_constexpr(tab_ltemp_max),
    tab_rho_min  = exp10_constexpr(tab_lrho_min),
    tab_rho_max  = exp10_constexpr(tab_lrho_max),
    tab_ltemp_delta = (tab_ltemp_max - tab_ltemp_min)/(double)(tab_ntemp - 1),
    tab_lrho_delta  = (tab_lrho_max - tab_lrho_min)/(double)(tab_nrho - 1),
    rho_extrapolation_margin = exp10_constexpr(tab_lrho_min
                             + tab_extrapolation_margin_irho*tab_lrho_delta),
    temp_extrapolation_margin = exp10_constexpr(tab_ltemp_min
                              + tab_extrapolation_margin_itemp*tab_ltemp_delta);

  static constexpr double
    eint_ele_deg_coef = exp10_constexpr(13.3369), // a constant in extrapolated electron degeneracy
    eint_ele_deg_thr1 = -7.50,          // threshold in ldin to switch to extrapolation
    eint_ele_deg_thr2 =  8.20,          // below this ltemp electron-positron contrib. is zero
    pres_ele_deg_coef = exp10_constexpr(12.4992), // extrapolated pressure for electron degenracy
    pres_ele_deg_coef2= exp10_constexpr(11.4198), // coefficient in a fit for el. degeneracy pressure
    // minimal value of the specific internal energy in the table (used in root finder)
    RGAS_CGS = 8.314462e7;              // gas constant (N_AVO * k_Boltzmann in cgs)

public:
  /**
  * @brief      Initialize tabulated EOS from Helmholtz
  *             Uses the path to EOS table (in ascii format).
  */
  static void init() {
    log_one(info) << "Reading tabulated EOS from file: "
                  << param::eos_tab_file_path << std::endl;
    eos_helm_init(param::eos_tab_file_path); // TODO: this does not scale
    log_one(info) << "... finished reading tabulated EoS." << std::endl;
  }

  /**
  * @brief      Compute pressure for tabulated EOS
  *               This function currently calculates pressure, soundspeed,
  *               energy density, entropy, etc. and sets the particle values
  *               in the function helm_eos_ptgiven. If you wish to expand stored
  *               values, add more particle.setValue(##) functions
  * @param      particle
  */
  static void compute_pressure(body & particle) {
    if (param::evolve_internal_energy) {
      compute_spct_given_rho_u(particle);
    } else {
      helm_eos_given_rho_t(particle);
    }
  } // compute_pressure_helm

  /**
  * @brief      Compute internal energy given temperature and density
  *
  * @param      particle
  */
  static void
  compute_internal_energy (body & particle) {
    const double rho = particle.getDensity(),
                temp = particle.getTemperature(),
                abar = particle.getAbar(),
                zbar = abar*particle.getElectronfraction();
    struct helm_eos_cache cache;
    helm_eos_update_cache(rho, abar, zbar, cache);
    double eint[5];
    get_eint_given_rho_temp(rho, temp, eint, cache);
    particle.setInternalenergy(eint[0]);
  } // compute_internal_energy


  /**
  * @brief      Compute speed of sound for tabulated EOS
  * @param      particle
  */
  static void
  compute_soundspeed(body & particle) {
    // Not used in current form
    const double pressure = particle.getPressure();
    const double entropy = particle.getEntropy();
    const double temp = particle.getTemperature();
    compute_spct_given_rho_u(particle);
    particle.setPressure(pressure);
    particle.setEntropy(entropy);
    particle.setTemperature(temp);
  } // compute_soundspeed

  /**
  * @brief      Compute pressure derivative of density for tabulated EOS
  * @param      particle
  */
  static double
  get_dPdrho(const body & particle) {
    //TODO: finish this function
    double dPdrho = 0;
    return dPdrho;
  }  // get_dPdrho

  static double get_dPdrhoInGeom(const body & particle) {
    return get_dPdrho(particle) / (C_LIGHT_CGS*C_LIGHT_CGS);
  }

  /**
  * @brief      Compute temperature
  *             Not used in current form
  *
  * @param      particle
  */
  static void compute_temperature(body & particle) {
    // TODO
  } // compute_temperature

  /**
  * @brief      Compute entropy, pressure, soundspeed and temperature
  *             given density and internal energy
  *
  * @param      particle
  * @returns    number of iterations (-1: root not bracketed)
  */
  static int
  root_finder_spct_given_rho_u (body & particle) {
    // particle data
    const double eint = particle.getInternalenergy()  // intergy: input
                      / phys::cfactor_from_<param::cgs_units>::energy
                      * phys::cfactor_from_<param::cgs_units>::mass,
                  rho = particle.getDensity()
                      / phys::cfactor_from_<param::cgs_units>::density,
                 abar = particle.getAbar(),
                 zbar = abar*particle.getElectronfraction();
    double temp = particle.getTemperature() // temperature: initial guess
                / phys::cfactor_from_<param::cgs_units>::temperature;
    const double tab_ltemp_factor = exp10(tab_ltemp_delta);
    int retval = 0; // return value: calls count to get_eint_given_rho_temp

    struct helm_eos_cache cache;
    helm_eos_update_cache(rho, abar, zbar, cache);
    const double
        KBOL = phys::kB
             / phys::cfactor_from_<param::cgs_units>::energy
             * phys::cfactor_from_<param::cgs_units>::temperature,
        AVO  = phys::NAvo,
        AR   = phys::arad
             / phys::cfactor_from_<param::cgs_units>::energy
             * phys::cfactor_from_<param::cgs_units>::volume
             * QU(phys::cfactor_from_<param::cgs_units>::temperature),
        ye   = zbar/abar,
        din  = cache.din,
        dxnidd = cache.dxnidd,
        ldin = cache.ldin;

    // arrays for the different contributions and their derivatives
    // all are initialized to zero
    double prad[5] = {0}, pion[5] = {0}, pele[5] = {0}, pcou[5] = {0};
    double erad[5] = {0}, eion[5] = {0}, eele[5] = {0}, ecou[5] = {0};
    double srad[5] = {0}, sion[5] = {0}, sele[5] = {0}, scou[5] = {0};
    double etot[5] = {0};
    // electron chemical potential and electron + positron number density
    double etaele[5] = {0}, xne[5] = {0};
    double p[5] = {0}, e[5] = {0}, s[5] = {0};

    if (cache.din < rho_extrapolation_margin
          || eint < tab_eint_min_at_rho(rho, ye)) {
        // handle low-T or low-rho extrapolation case
        double edeg = ye*((ldin < eint_ele_deg_thr1)
                    ? eint_ele_deg_coef*rho
                    : exp10(eint_ele_deg(ldin)));
        double temp_max_estimate = sqrt(sqrt(rho*(edeg - eint)/AR));
        double temp1 = temp_max_estimate, temp2;
        int niter = 0;
        do {
            temp2 = temp1;
            get_eint_given_rho_temp(rho, temp2, etot, cache);
            temp1 = temp2 - (etot[VALUE] - eint)/etot[DTEMP];
            temp1 = std::max(temp1, temp_extrapolation_margin);
        } while (std::abs(1 - temp1/temp2) > HELM_EOS_EPS && niter++ < 50);
        //helm_eos_rad(rho, temp1, prad, erad, srad);
        //helm_eos_ion(rho, temp1, pion, eion, sion, cache);
        //helm_eos_ele_offtab(rho, temp1, pele, eele, sele, etaele, xne, cache);
        //printf("%24.15e %24.15e %24.15e\n", rho, etot[0], erad[0] + eion[0] + eele[0]);
        temp = temp1;
        helm_eos_rad(rho, temp, prad, erad, srad);
        helm_eos_ion(rho, temp, pion, eion, sion, cache);
        helm_eos_ele(rho, temp, pele, eele, sele, etaele, xne, cache);
        helm_eos_cou(rho, temp, pcou, ecou, scou, cache);

        p[VALUE] = prad[VALUE] + pion[VALUE] + pele[VALUE] + pcou[VALUE];
        e[VALUE] = erad[VALUE] + eion[VALUE] + eele[VALUE] + ecou[VALUE];
        s[VALUE] = srad[VALUE] + sion[VALUE] + sele[VALUE] + scou[VALUE];
        p[DTEMP] = prad[DTEMP] + pion[DTEMP] + pele[DTEMP] + pcou[DTEMP];
        p[DRHO]  = prad[DRHO]  + pion[DRHO]  + pele[DRHO]  + pcou[DRHO];
        e[DTEMP] = erad[DTEMP] + eion[DTEMP] + eele[DTEMP] + ecou[DTEMP];

        const double C_LIGHT_CGS = phys::clight
                                 / phys::cfactor_from_<param::cgs_units>::velocity;

        double cv       = e[DTEMP];
        double chit     = temp/p[VALUE]*p[DTEMP];
        double chid     = p[DRHO]*rho/p[VALUE];
        double x        = p[VALUE]/rho * chit/(temp*cv);
        double gamma_1  = chit*x + chid;
        double sound    = C_LIGHT_CGS*sqrt(gamma_1
                        /(1 + (e[VALUE] + SQ(C_LIGHT_CGS))*rho/p[VALUE]));

        particle.setSoundspeed(
                 sound * phys::cfactor_from_<param::cgs_units>::velocity);
        particle.setPressure(
                 p[VALUE] * phys::cfactor_from_<param::cgs_units>::pressure);
        particle.setEntropy(s[VALUE]); //TODO: you forgot to add conversion factor for the entropy!
        particle.setTemperature(
                 temp * phys::cfactor_from_<param::cgs_units>::temperature);

        return niter;
    }

    // initial guess: try to bracket the temperature from previous values
    int jat = (log10(temp) - tab_ltemp_min)/tab_ltemp_delta;

    double temp1 = tab_temp_min*exp10(jat*tab_ltemp_delta);
    get_eint_given_rho_temp(rho, temp1, etot, cache);
    retval++;
    double e1 = etot[VALUE];
    double dedt1 = etot[DTEMP];

    double temp2 = tab_temp_min*exp10((jat + 1)*tab_ltemp_delta);
    get_eint_given_rho_temp(rho, temp2, etot, cache);
    retval++;
    double e2 = etot[VALUE];
    double dedt2 = etot[DTEMP];

    if (eint < e1 && jat > 0) {
      // test the neighboring grid cell below
      --jat;
      temp2 = temp1;
      e2 = e1;
      dedt2 = dedt1;

      temp1 = tab_temp_min*exp10(jat*tab_ltemp_delta);
      get_eint_given_rho_temp(rho, temp1, etot, cache);
      retval++;
      e1 = etot[VALUE];
      dedt1 = etot[DTEMP];
    }

    if (eint > e2 && jat < tab_ntemp - 1) {
      // test the neighboring grid cell above
      ++jat;
      temp1 = temp2;
      e1 = e2;
      dedt1 = dedt2;

      temp2 = tab_temp_min*exp10((jat + 1)*tab_ltemp_delta);
      get_eint_given_rho_temp(rho, temp2, etot, cache);
      retval++;
      e2 = etot[VALUE];
      dedt2 = etot[DTEMP];
    }

    // bisection on the temperature grid
    int j1 = jat, j2 = jat + 1;
    if ((eint - e1)*(eint - e2) > 0.) {

      if (eint < e1) {
        j2 = j1;
        temp2 = temp1;
        dedt2 = dedt1;

        j1 = 0;
        temp1 = tab_temp_min;
        get_eint_given_rho_temp(rho, temp1, etot, cache);
        retval++;
        e1 = etot[VALUE];
        dedt1 = etot[DTEMP];
      }

      if (eint > e2) {
        j1 = j2;
        temp1 = temp2;
        dedt1 = dedt2;

        j2 = tab_ntemp - 1;
        temp2 = tab_temp_max;
        get_eint_given_rho_temp(rho, temp2, etot, cache);
        retval++;
        e2 = etot[VALUE];
        dedt2 = etot[DTEMP];
      }

      if ((e1 - eint)*(e2 - eint) > 0) {
        // root not bracketed
        return -1;
      }

      // bisection
      const int itmax = floor(log(tab_ntemp)/log(2.)) + 1;
      double temp12 = temp1;
      int it;
      for (it = 0; it < itmax; ++it) {
        if (j2 - j1 < 2)
          break;
        jat = (j1 + j2)/2;
        temp12 = tab_temp_min*exp10(jat*tab_ltemp_delta);
        get_eint_given_rho_temp(rho, temp12, etot, cache);
        retval++;
        double e12 = etot[VALUE];
        double dedt12 = etot[DTEMP];
        if (eint < e12) {
          temp2 = temp12;
          e2 = e12;
          dedt2 = dedt12;
          j2 = jat;
        }
        else {
          temp1 = temp12;
          e1 = e12;
          dedt1 = dedt12;
          j1 = jat;
        }
      } // it = 0..itmax
    } // if not (e1 < eint < e2)

    // At this point, we have found a single grid cell [jat, jat+1]
    // that brackets the temperature: temp \in [temp1, temp2],
    // where temp1 = temp_grid[jat], temp2 = temp_grid[jat+1].

    assert ((eint - e1)*(eint - e2) <= 0.);
    temp = temp1;

    // check if either of end points is the root
    if (std::abs(eint - e1) < std::abs(e1)*HELM_EOS_EPS) {
      // do nothing: temp is already == temp1
    }
    else if (std::abs(eint - e2) < std::abs(e2)*HELM_EOS_EPS) {
      temp = temp2;
    }
    else {
      for (int nr = 0; nr < HELM_EOS_MAXITER; ++nr) {
        double temp_p = temp; // temperature from previous iteration
        if ((e2 - e1)*(dedt2 - dedt1) > 0)
          // use tangent at the right endpoint
          temp = temp2 - (e2 - eint)/dedt2;
        else
          // use tangent at the left endpoint
          temp = temp1 - (e1 - eint)/dedt1;

        if (temp <= temp1 or temp >= temp2) {
          // if root outside interval, use secant
          temp = temp1 + (eint - e1)*(temp2 - temp1)/(e2 - e1);

          // secant can be quite slow;
          // check if relative change in temp is too tiny and if so,
          // use bisection
          if (std::abs(temp - temp_p)/(temp2 - temp1) < 0.1)
            temp = (temp1 + temp2)*.5;
        }

        get_eint_given_rho_temp(rho, temp, etot, cache);
        retval++;
        if (std::abs(etot[VALUE] - eint) < std::abs(eint)*HELM_EOS_EPS)
          break;

        // bisection steps
        if (etot[VALUE] < eint) {
          e1 = etot[VALUE];
          temp1 = temp;
          dedt1 = etot[DTEMP];
        }
        else {
          e2 = etot[VALUE];
          temp2 = temp;
          dedt2 = etot[DTEMP];
        }

      } // nr
    }

    helm_eos_rad(rho, temp, prad, erad, srad);
    helm_eos_ion(rho, temp, pion, eion, sion, cache);
    helm_eos_ele(rho, temp, pele, eele, sele, etaele, xne, cache);
    helm_eos_cou(rho, temp, pcou, ecou, scou, cache);

    p[VALUE] = prad[VALUE] + pion[VALUE] + pele[VALUE] + pcou[VALUE];
    e[VALUE] = erad[VALUE] + eion[VALUE] + eele[VALUE] + ecou[VALUE];
    s[VALUE] = srad[VALUE] + sion[VALUE] + sele[VALUE] + scou[VALUE];
    p[DTEMP] = prad[DTEMP] + pion[DTEMP] + pele[DTEMP] + pcou[DTEMP];
    p[DRHO]  = prad[DRHO]  + pion[DRHO]  + pele[DRHO]  + pcou[DRHO];
    e[DTEMP] = erad[DTEMP] + eion[DTEMP] + eele[DTEMP] + ecou[DTEMP];

    const double C_LIGHT_CGS = phys::clight
                             / phys::cfactor_from_<param::cgs_units>::velocity;

    double cv       = e[DTEMP];
    double chit     = temp/p[VALUE]*p[DTEMP];
    double chid     = p[DRHO]*rho/p[VALUE];
    double x        = p[VALUE]/rho * chit/(temp*cv);
    double gamma_1  = chit*x + chid;
    double sound    = C_LIGHT_CGS*sqrt(gamma_1
                    /(1 + (e[VALUE] + SQ(C_LIGHT_CGS))*rho/p[VALUE]));

    particle.setSoundspeed(
             sound * phys::cfactor_from_<param::cgs_units>::velocity);
    particle.setPressure(
             p[VALUE] * phys::cfactor_from_<param::cgs_units>::pressure);
    particle.setEntropy(s[VALUE]); //TODO: you forgot to add conversion factor for the entropy!
    particle.setTemperature(
             temp * phys::cfactor_from_<param::cgs_units>::temperature);

//    if (cache.din < rho_extrapolation_margin) {
//     printf("%24.15e: %24.15e\n", temp, (1. - temp_max_estimate/temp));
//    }
//helm_eos_ele(rho, temp, pele, eele, sele, etaele, xne, cache);
//printf("%15.9e  %15.9e  %24.15e\n", rho, temp, eele[VALUE]);

    return retval;

  } // root_finder_spct_given_rho_u


  /**
  * @brief      Wrapper for the above with void return type
  *
  * @param      particle
  * @returns    number of iterations (-1: root not bracketed)
  */
  static void
  compute_spct_given_rho_u (body & particle) {
    int ncounts = root_finder_spct_given_rho_u (particle);
    if (ncounts < 0) {
      log_one(warn) << "failed to invert at {rho, temp} = "
        << std::scientific << std::setprecision(20)
        << "{" << particle.getDensity()
        << ", " << particle.getTemperature()
        << "} from internal energy "
        << particle.getInternalenergy() << std::endl;
    }
  }


  /////////////////////////////////////////////////////////////////////////////
  // GETTING ENTROPY (S) FROM INITIAL CONDITIONS GIVEN RHO AND EINT
  static void
  compute_entropy(body & particle) {
    const double pressure = particle.getPressure();
    const double soundspeed = particle.getSoundspeed();
    const double temp = particle.getTemperature();
    compute_spct_given_rho_u(particle);
    particle.setPressure(pressure);
    particle.setSoundspeed(soundspeed);
    particle.setTemperature(temp);
  } // compute_entropy

  /////////////////////////////////////////////////////////////////////////////
  // GETTING ENTROPY (S) FROM INITIAL CONDITIONS GIVEN TEMP AND RHO
  static void
  compute_entropy_given_rho_temp(body & particle) {
    double temp = particle.getTemperature();
    const double rho = particle.getDensity(),
                abar = particle.getAbar(),
                zbar = abar*particle.getElectronfraction();
    struct helm_eos_cache cache;

    helm_eos_update_cache(rho, abar, zbar, cache);

    temp = std::min(tab_temp_max, std::max(tab_temp_min, temp));
    double ent = entropy_helm_eos_rad(rho, temp)
               + entropy_helm_eos_ion(rho, temp, cache);
               + entropy_helm_eos_ele(rho, temp, cache);
               + entropy_helm_eos_cou(rho, temp, cache);
    particle.setEntropy(ent);
  } // compute_entropy_given_rho_temp

  /**
  * @brief      Perform consistency checks on the equation of state
  */
  static void
  root_finder_check() {

    // begin table solve
    struct helm_eos_cache cache;
    double abar = 13.;
    double zbar = 6.;
    body particle;
    log_one(info) << "Helmholtz EoS consistency check" << std::endl;
    //printf ("# 1:i 2:j 3:rho 4:temp 5:eint 6:entropy 7:pressure\n");
    double eint_sum = 0.0;

    particle.setAbar(abar);
    particle.setElectronfraction(zbar/abar);

    double prad[5] = {0}, pion[5] = {0}, pele[5] = {0}, pcou[5] = {0};
    double erad[5] = {0}, eion[5] = {0}, eele[5] = {0}, ecou[5] = {0};
    double srad[5] = {0}, sion[5] = {0}, sele[5] = {0}, scou[5] = {0};
    double etaele[5] = {0}, xne[5] = {0};

    /*
    const int
      test_nrho = 1,
      test_ntemp = 1000,
      ninv = test_nrho*test_ntemp;

    const double
      test_lrho_min = 0.,//  + log10(abar/zbar),
      test_lrho_max = test_lrho_min,
      test_ltemp_min = 1.,
      test_ltemp_max = 12.,
      test_lrho_delta = (test_lrho_max - test_lrho_min)
                      / std::max(1,test_nrho - 1),
      test_ltemp_delta = (test_ltemp_max - test_ltemp_min)
                      / std::max(1,test_ntemp - 1);

    const int
      test_nrho = 1000,
      test_ntemp = 1,
      ninv = test_nrho*test_ntemp;

    const double
      test_lrho_min = -22. + log10(abar/zbar),
      test_lrho_max =  15. + log10(zbar/abar),
      test_ltemp_min = log10(273.15),
      test_ltemp_max = test_ltemp_min,
      test_lrho_delta = (test_lrho_max - test_lrho_min)
                      / std::max(1,test_nrho - 1),
      test_ltemp_delta = (test_ltemp_max - test_ltemp_min)
                      / std::max(1,test_ntemp - 1);

    //// NEGATIVE SLOPE CASE:
    // abar = 48.; zbar = 23.;
    // rho = 2.71227257933202126878e+04;
    // double temp = 3e+7;
    */

    const int
      test_nrho = 1503,
      test_ntemp = 1710,
      ninv = test_nrho*test_ntemp;

    const double
      test_lrho_min = -22.0 + log10(abar/zbar),
      test_lrho_max =  7.0 + log10(zbar/abar),
      test_ltemp_min = -1.,
      test_ltemp_max = 11.999,
      test_lrho_delta = (test_lrho_max - test_lrho_min)/(test_nrho - 1),
      test_ltemp_delta = (test_ltemp_max - test_ltemp_min)/(test_ntemp - 1);

    double temp_guess = 79999.;
    log_one(info) << std::endl
      << "Testing {rho, T} <--> {rho, eint} inversions"
      << " on a log-uniform grid of values (" << ninv << " inversions)"
      << std::endl
      << "Parameters: " << std::endl
      << " - Abar = " << abar << std::endl
      << " - Zbar = " << zbar << std::endl
      << " - tolerance (HELM_EOS_EPS) = " << HELM_EOS_EPS << std::endl
      << " - max. interations (HELM_EOS_MAXITER) = " << HELM_EOS_MAXITER
      << std::endl
      << " - log10(rho)  in {"<<test_lrho_min << ", " << test_lrho_max<<"}" << std::endl
      << " - log10(temp) in {"<<test_ltemp_min<<", "<<test_ltemp_max<<"}"   << std::endl
      << " - grid: {N_rho x N_temp} = {"<<test_nrho<<" x "<<test_ntemp<<"}" << std::endl
      << " - temperature guess: " << temp_guess
      << std::endl;
    double pres_L2_error = 0., pres_Lmax_error = 0.;
    double temp_L2_error = 0., temp_Lmax_error = 0.;
    double entr_L2_error = 0., entr_Lmax_error = 0.;
    double ncalls_ave = 0.;
    double pres_max_error_rho = 0., pres_max_error_temp = 0.;
    double temp_max_error_rho = 0., temp_max_error_temp = 0.;
    double entr_max_error_rho = 0., entr_max_error_temp = 0.;
    int successful_inversions_count = 0;
    for (int i = 0; i < test_nrho; i++) {
      double rho = exp10(test_lrho_min + i*test_lrho_delta);
      helm_eos_update_cache(rho, abar, zbar, cache);
      for (int j = 0; j < test_ntemp; j++) {
        double temp = exp10(test_ltemp_min + j*test_ltemp_delta);

        helm_eos_rad(rho, temp, prad, erad, srad);
        helm_eos_ion(rho, temp, pion, eion, sion, cache);
        helm_eos_ele_offtab(rho, temp, pele, eele, sele, etaele, xne, cache);
        double eele_0 = eele[0];
        double eele_1 = eele[1];
        double eele_2 = eele[2];
        double pele_0 = pele[0];
        double pele_1 = pele[1];
        double pele_2 = pele[2];
        helm_eos_ele(rho, temp, pele, eele, sele, etaele, xne, cache); //, false);
        helm_eos_cou(rho, temp, pcou, ecou, scou, cache);
        double eint = erad[0] + eion[0] + eele[0] + ecou[0];
        double entr = srad[0] + sion[0] + sele[0] + scou[0];
        double pres = prad[0] + pion[0] + pele[0] + pcou[0];

        double rho_units
               = rho * phys::cfactor_from_<param::cgs_units>::density;
        double eint_units
               = eint
               * phys::cfactor_from_<param::cgs_units>::energy
               / phys::cfactor_from_<param::cgs_units>::mass;
        double temp_guess_units
               = temp_guess
               * phys::cfactor_from_<param::cgs_units>::temperature;

        particle.setDensity(rho_units);
        particle.setInternalenergy(eint_units);
        particle.setTemperature(temp_guess_units); // initial guess
        int ncalls = root_finder_spct_given_rho_u (particle);

        if (ncalls < 0) {
          log_one(warn) << "failed to invert at {rho, temp} = "
            << std::scientific << std::setprecision(20)
            << "{" << rho << ", " << temp << "}" << std::endl;
        }
        else {
          successful_inversions_count++;
          ncalls_ave += (double)ncalls;

          double temp_units
                 = particle.getTemperature()
                 / phys::cfactor_from_<param::cgs_units>::temperature;
          double eps = 1.0 - temp/temp_units;
          temp_L2_error += eps*eps;
          if (temp_Lmax_error < std::abs(eps)) {
            temp_Lmax_error = std::abs(eps);
            temp_max_error_temp = temp;
            temp_max_error_rho  = rho;
          }

          double pres_units
                 = particle.getPressure()
                 / phys::cfactor_from_<param::cgs_units>::pressure;
          eps = 1.0 - pres/pres_units;
          pres_L2_error += eps*eps;
          if (pres_Lmax_error < std::abs(eps)) {
            pres_Lmax_error = std::abs(eps);
            pres_max_error_temp = temp;
            pres_max_error_rho  = rho;
          }

          eps = 1.0 - entr/particle.getEntropy();
          entr_L2_error += eps*eps;
          if (entr_Lmax_error < std::abs(eps)) {
            entr_Lmax_error = std::abs(eps);
            entr_max_error_temp = temp;
            entr_max_error_rho  = rho;
          }
        }

//printf ("%14.7e  %14.7e   %24.17e  %24.17e  %24.17e    %24.17e  %24.17e  %24.17e\n",
//         rho,temp, eele[0],eele[1],eele[2],
//         eele_0, eele_1, eele_2);
//         rho,temp, eele[0],eele[1],eele[2], eele_0, eele_1, eele_2);
//printf ("%14.7e  %14.7e   %24.17e  %24.17e  %24.17e    %24.17e  %24.17e  %24.17e\n",
//         rho,temp, pele[0],pele[1],pele[2], pele_0, pele_1, pele_2);

      } // j: temperature index

      //// empty line for gnuplot-friendly output
      //std::cout << std::endl;
    } // i: density index

    temp_L2_error = sqrt(temp_L2_error
                  / std::max(1, successful_inversions_count));
    pres_L2_error = sqrt(pres_L2_error
                  / std::max(1, successful_inversions_count));
    entr_L2_error = sqrt(entr_L2_error
                  / std::max(1, successful_inversions_count));
    ncalls_ave /= (double)std::max(1, successful_inversions_count);

    log_one(info) << std::endl
      << "Summary:" << std::endl
      << " - successful inversions: " << successful_inversions_count
      << " out of " << ninv << ";" << std::endl
      << " - average number of iterations per inversion: "
      << std::fixed << std::setprecision(1)
      << ncalls_ave << std::endl
      << std::scientific << std::setprecision(4)
      << " - temperature: average rms error = " << temp_L2_error
      << ", max error = " << temp_Lmax_error
      << std::scientific << std::setprecision(20)
      << " at {rho, temp} = {" << temp_max_error_rho
      << ", " << temp_max_error_temp << "}" << std::endl
      << std::scientific << std::setprecision(4)
      << " - pressure:    average rms error = " << pres_L2_error
      << ", max error = " << pres_Lmax_error
      << std::scientific << std::setprecision(20)
      << " at {rho, temp} = {" << pres_max_error_rho
      << ", " << pres_max_error_temp << "}" << std::endl
      << std::scientific << std::setprecision(4)
      << " - entropy:     average rms error = " << entr_L2_error
      << ", max error = " << entr_Lmax_error
      << std::scientific << std::setprecision(20)
      << " at {rho, temp} = {" << entr_max_error_rho
      << ", " << entr_max_error_temp << "}" << std::endl;

  } // consistency_check

  /**
  * @brief      Perform various checks on the table
  */
  static void
  table_check() {

    // begin table solve
    struct helm_eos_cache cache;
    double abar = 1.;
    double zbar = 1.;
    body particle;
    printf ("# Helmholtz EoS table check\n");
    printf ("# 1:i 2:j 3:rho 4:temp 5:F\n");

    particle.setAbar(abar);
    particle.setElectronfraction(zbar/abar);

    double prad[5] = {0}, pion[5] = {0}, pele[5] = {0}, pcou[5] = {0};
    double erad[5] = {0}, eion[5] = {0}, eele[5] = {0}, ecou[5] = {0};
    double srad[5] = {0}, sion[5] = {0}, sele[5] = {0}, scou[5] = {0};
    double etaele[5] = {0}, xne[5] = {0};

    double dftemp = exp10(tab_ltemp_delta);

    for (int i = 0; i < tab_nrho; i++) {
      double rho = tab_rho_min*exp10(i*tab_lrho_delta);
      for (int j = 0; j < tab_ntemp; j++) {
        double temp = tab_temp_min*exp10(j*tab_ltemp_delta);

        double free_en = helm_eos_table_ptr->f[i][j];
        double df_t = helm_eos_table_ptr->ft[i][j];
        double eint = free_en - temp*df_t;
        double dedt = -temp*helm_eos_table_ptr->ftt[i][j];

        printf ("% 3d  % 3d  %24.17e  %24.17e  %24.17e  %24.17e\n",
                   i,    j,   rho,     temp,    eint,   dedt);

        if (i > 0) {
          double free_en1 = helm_eos_table_ptr->f[i-1][j];
          double df_t1 = helm_eos_table_ptr->ft[i-1][j];
          double eint1 = free_en - temp*df_t;

          if (eint < eint1) {
            log_one(error) << "internal energy decreasing at (i-1,j) = ("
                           << (i-1) << ", " << j << "): "
                           << eint1 << " --> " << eint << std::endl;
          }
        }

        if (j > 0) {
          double free_en1 = helm_eos_table_ptr->f[i][j-1];
          double df_t1 = helm_eos_table_ptr->ft[i][j-1];
          double eint1 = free_en - temp/dftemp*df_t;

          if (eint < eint1) {
            log_one(error) << "internal energy decreasing at (i,j-1) = ("
                           << i << ", " << (j-1) << "): "
                           << std::scientific << std::setprecision(16)
                           << eint1 << " --> " << eint << std::endl;
          }
        }

        } // j: temperature index
      // empty line for gnuplot output
      std::cout << std::endl;
    } // i: density index

  } // table_check

  // TODO
  static get_quantity_t get_dpdrho_at_temp;

private:
  /**
   * @brief      creates structure for equation of state quantities and its derivatives
   */
  struct state_value {
    double val;   // value of eos quantity
    double drho;  // derivative of quantity w.r.t. density
    double dtemp; // derivative of quantity w.r.t. temperature
    double dabar; // derivative of quantity w.r.t. abar
    double dzbar; // derivative of quantity w.r.t. zbar
  };

  /**
   * @brief      creates structure for results from the eos calculation
   */
  struct eos_result {
    double temp;       // value of temperature
    double p[5];       // pressure properties
    double e[5];       // specific energy properties
    double s[5];       // specific entropy properties
    double etaele[5];  // degeneracy parameter: electron chemical potential
    double nep[5];     // electro+positron number density properties
    double cv;         // specific heat at constant volume
    double cp;         // specific heat at constant pressure
    double chit;       // temperature exponent from Cox & Giuli
    double chid;       // density exponent from Cox & Giuli
    double gamma_1, gamma_2, gamma_3; // gamma parameters from Cox & Giuli
    double nabla_ad;   // adiabatic nabla
    double delta, phi; // parameters from Kippenhahn & Weigert (6.6)
    double sound;      // relativistic sound speed
    double abar;       // mean mass number
    double zbar;       // mean charge number
  };

  /**
   * @brief      creates cache to store pertinent values
   */
  struct helm_eos_cache {
    double abar, zbar; // mean mass and charge number
    double ye;         // zbar/abar
    double din, ldin;  // rho*ye and log
    double ytot;       // 1.0/abar;
    double xni;        // AVO*ytot*rho
    double dxnidd;     // AVO*ytot
    double dxnida;     // -xni*ytot
  };

  /**
   * @brief      declare structures TODO: move to eos_utils?
   */

  typedef double helm_eos_table_entry[tab_nrho][tab_ntemp];

  static struct interpolating_function_1d
    eint_ele_deg,  // degeneracy part
    dedd_ele_deg,  // degeneracy part: derivative wrt density
    eint_ele_ep,   // electron-positron pairs
    dedt_ele_ep,   // electron-positron pairs: derivative wrt temperature
    pres_ele_deg,  // pressure: degeneracy part
    dpdd_ele_deg,  // pressure: degeneracy part, derivative wrt density
    pres_ele_ep,   // pressure: electron-positron pairs
    dpdt_ele_ep;   // pressure: electron-positron pairs, derivative wrt temp

  /**
   * @brief      creates structure for storing the helmholtz datafile
   */
  struct helm_eos_table {
    // density and temperature ranges
    double temp[tab_ntemp];
    double rho[tab_nrho];
    // d means derivative w.r.t. density; t means derivative w.r.t. temperature
    // Helmholtz free energy
    helm_eos_table_entry f;
    helm_eos_table_entry fd;
    helm_eos_table_entry ft;
    helm_eos_table_entry fdd;
    helm_eos_table_entry ftt;
    helm_eos_table_entry fdt;
    helm_eos_table_entry fddt;
    helm_eos_table_entry fdtt;
    helm_eos_table_entry fddtt;
    // pressure derivative w.r.t. density
    helm_eos_table_entry dpdf;
    helm_eos_table_entry dpdfd;
    helm_eos_table_entry dpdft;
    helm_eos_table_entry dpdfdd;
    helm_eos_table_entry dpdftt;
    helm_eos_table_entry dpdfdt;
    // chemical potential
    helm_eos_table_entry ef;
    helm_eos_table_entry efd;
    helm_eos_table_entry eft;
    helm_eos_table_entry efdd;
    helm_eos_table_entry eftt;
    helm_eos_table_entry efdt;
    // number density
    helm_eos_table_entry xf;
    helm_eos_table_entry xfd;
    helm_eos_table_entry xft;
    helm_eos_table_entry xfdd;
    helm_eos_table_entry xftt;
    helm_eos_table_entry xfdt;
    // species information
    //int nspecies;
    //double *na;
    //double *nai;
    //double *nz;
  };

  static helm_eos_table* helm_eos_table_ptr;

  /////////////////////////////////////////////////////////////////////////////
  // QUINTIC HERMITE POLYNOMIALS
  /////////////////////////////////////////////////////////////////////////////
  // PSI0 AND ITS DERIVATIVES
  static double psi0(const double z) {
    return CU(z) * (z * (-6.0 * z + 15.0) - 10.0) + 1.0;
  }
  static double dpsi0(const double z) {
    return SQ(z) * (z * (-30.0 * z + 60.0) - 30.0);
  }
  static double ddpsi0(const double z) {
    return z * (z * (-120.0 * z + 180.0) - 60.0);
  }

  /////////////////////////////////////////////////////////////////////////////
  // PSI1 AND ITS DERIVATIVES
  static double psi1(const double z) {
    return z * (SQ(z) * (z * (-3.0 * z + 8.0) - 6.0) + 1.0);
  }
  static double dpsi1(const double z) {
    return SQ(z) * (z * (-15.0 * z + 32.0) - 18.0) + 1.0;
  }
  static double ddpsi1(const double z) {
    return z * (z * (-60.0 * z + 96.0) - 36.0);
  }

  /////////////////////////////////////////////////////////////////////////////
  // PSI2 AND ITS DERIVATIVES
  static double psi2(const double z) {
    return 0.5 * SQ(z) * (z * (z * (-z + 3.0) - 3.0) + 1.0);
  }
  static double dpsi2(const double z) {
    return 0.5 * z * ( z * (z * (-5.0 * z + 12.0) - 9.0) + 2.0);
  }
  static double ddpsi2(const double z) {
    return 0.5 * (z * (z * (-20.0 * z + 36.0) - 18.0) + 2.0);
  }

  /////////////////////////////////////////////////////////////////////////////
  // BIQUINTIC HERMITE POLYNOMIAL
  static double h5xy(const double fi[36],
    const double w0t,   const double w1t,    const double w2t,
    const double w0mt,  const double w1mt,   const double w2mt,
    const double w0d,   const double w1d,    const double w2d,
    const double w0md,  const double w1md,   const double w2md) {
    return
      fi[0] *w0d*w0t + fi[1] *w0md*w0t + fi[2] *w0d*w0mt + fi[3] *w0md*w0mt
    + fi[4] *w0d*w1t + fi[5] *w0md*w1t + fi[6] *w0d*w1mt + fi[7] *w0md*w1mt
    + fi[8] *w0d*w2t + fi[9] *w0md*w2t + fi[10]*w0d*w2mt + fi[11]*w0md*w2mt
    + fi[12]*w1d*w0t + fi[13]*w1md*w0t + fi[14]*w1d*w0mt + fi[15]*w1md*w0mt
    + fi[16]*w2d*w0t + fi[17]*w2md*w0t + fi[18]*w2d*w0mt + fi[19]*w2md*w0mt
    + fi[20]*w1d*w1t + fi[21]*w1md*w1t + fi[22]*w1d*w1mt + fi[23]*w1md*w1mt
    + fi[24]*w2d*w1t + fi[25]*w2md*w1t + fi[26]*w2d*w1mt + fi[27]*w2md*w1mt
    + fi[28]*w1d*w2t + fi[29]*w1md*w2t + fi[30]*w1d*w2mt + fi[31]*w1md*w2mt
    + fi[32]*w2d*w2t + fi[33]*w2md*w2t + fi[34]*w2d*w2mt + fi[35]*w2md*w2mt;
  }

  /////////////////////////////////////////////////////////////////////////////
  // CUBIC HERMITE POLYNOMIALS
  /////////////////////////////////////////////////////////////////////////////
  // PSI0 AND ITS DERIVATIVE
  static double xpsi0(double z) {
    return SQ(z) * (2.0 * z - 3.0) + 1.0;
  }
  static double xdpsi0(double z) {
    return z * (6.0 * z - 6.0);
  }

  /////////////////////////////////////////////////////////////////////////////
  // PSI1 AND ITS DERIVATIVE
  static double xpsi1(double z) {
    return z * (z * (z - 2.0) + 1.0);
  }
  static double xdpsi1(double z) {
    return z * (3.0 * z - 4.0) + 1.0;
  }

  /////////////////////////////////////////////////////////////////////////////
  // CUBIC HERMITE POLYNOMIAL FOR ONE VARIABLE
  static double h3x(const double x, const double x1, const double dx,
      const double fi[4]) {

    const double
      z = (x - x1)/dx,
      z1 = 1. - z;

    return fi[0]*xpsi0(z) + fi[1]*xpsi0(z1)
     + dx*(fi[2]*xpsi1(z) - fi[3]*xpsi1(z1));
  }

  /////////////////////////////////////////////////////////////////////////////
  // FIRST DERIVATIVE OF THE CUBIC HERMITE POLYNOMIAL FOR ONE VARIABLE
  static double d_h3x(const double x, const double x1, const double dx,
      const double fi[4]) {

    const double
      z = (x - x1)/dx,
      z1 = 1. - z;

    return dx*(fi[0]*xdpsi0(z) + fi[1]*xdpsi0(z1)
         + dx*(fi[2]*xdpsi1(z) - fi[3]*xdpsi1(z1)));
  }

  /////////////////////////////////////////////////////////////////////////////
  // BICUBIC HERMITE POLYNOMIAL
  static double h3xy(const double fi[16],
      const double w0t, const double w1t, const double w0mt, const double w1mt,
      const double w0d, const double w1d, const double w0md, const double w1md) {
    return
      fi[0] *w0d*w0t + fi[1] *w0md*w0t + fi[2] *w0d*w0mt + fi[3] *w0md*w0mt
    + fi[4] *w0d*w1t + fi[5] *w0md*w1t + fi[6] *w0d*w1mt + fi[7] *w0md*w1mt
    + fi[8] *w1d*w0t + fi[9] *w1md*w0t + fi[10]*w1d*w0mt + fi[11]*w1md*w0mt
    + fi[12]*w1d*w1t + fi[13]*w1md*w1t + fi[14]*w1d*w1mt + fi[15]*w1md*w1mt;
  }

  /////////////////////////////////////////////////////////////////////////////
  // HERE ARE THE FUNCTIONS THAT CALCULATE THE DIFFERENT CONTRIBUTIONS OF THE EOS
  //   THE FIVE-ELEMENT ARRAYS CORRESPOND TO THE ACTUAL QUANTITY and
  //   ITS DERIVATIVES W.R.T. DENSITY, TEMPERATURE, ABAR, AND ZBAR
  /////////////////////////////////////////////////////////////////////////////
  // RADIATION SECTION
  static void
  helm_eos_rad(const double rho, const double temp,
      double prad[5], double erad[5], double srad[5]) {
    const double rhoi  = 1.0 / rho;
    const double tempi = 1.0 / temp;
    const double AR = phys::arad
                    / phys::cfactor_from_<param::cgs_units>::energy
                    * phys::cfactor_from_<param::cgs_units>::volume
                    * QU(phys::cfactor_from_<param::cgs_units>::temperature);

    prad[0] = AR / 3. * QU(temp);           // prad
    prad[1] = 0.;                           // dprad dd
    prad[2] = 4. * prad[0] * tempi;         // dprad dt
    prad[3] = 0.;                           // dprad da
    prad[4] = 0.;                           // dprad dz

    erad[0] = 3. * prad[0] * rhoi;          // erad
    erad[1] =-erad[0] * rhoi;               // derad dd
    erad[2] = 3. * prad[2] * rhoi;          // derad dt
    erad[3] = 0.;                           // derad da
    erad[4] = 0.;                           // derad dz

    srad[0] = 4. * prad[0] * rhoi * tempi;  // srad
    srad[1] =-srad[0] * rhoi;               // dsrad dd
    srad[2] = 4. * AR * rhoi * temp * temp; // dsrad dt
    srad[3] = 0.;                           // dsrad da
    srad[4] = 0.;                           // dsrad dz
  } //helm_eos_rad

  /////////////////////////////////////////////////////////////////////////////
  // ION SECTION
  static void
  helm_eos_ion(const double rho, const double temp,
      double pion[5], double eion[5], double sion[5],
      const struct helm_eos_cache & cache) {
    const double KBOL = phys::kB
                      / phys::cfactor_from_<param::cgs_units>::energy
                      * phys::cfactor_from_<param::cgs_units>::temperature;
    const double AMU  = phys::amu
                      / phys::cfactor_from_<param::cgs_units>::mass;
    const double HPL  = phys::hplanck
                      / phys::cfactor_from_<param::cgs_units>::energy
                      / phys::cfactor_from_<param::cgs_units>::time;
    const double AVO = phys::NAvo;
    const double kt  = KBOL * temp;
    const double abar = cache.abar;
    const double ytot = cache.ytot;
    const double xni  = cache.xni;
    const double dxnidd = cache.dxnidd;
    const double dxnida = cache.dxnida;
    //const double y = cache.ywot + cache.lswot15 + 1.5 * ltemp;
    const double s = (2.0 * M_PI * AMU * kt) / SQ(HPL);
    const double y = log(abar * abar * sqrt(abar) / (rho * AVO) * s * sqrt(s));

    pion[0] = xni * kt;                   // pion
    pion[1] = dxnidd * kt;                // dpion dd
    pion[2] = xni * KBOL;                 // dpion dt
    pion[3] = dxnida * kt;                // dpion da
    pion[4] = 0.0;                        // dpion dz

    eion[0] = 1.5 * pion[1];              // eion
    eion[1] = 0.0;                        // deion dd
    eion[2] = 1.5 * dxnidd * KBOL;        // deion dt
    eion[3] =-eion[0] * ytot;             // deion da
    eion[4] = 0.0;                        // deion dz

    sion[0] = KBOL * dxnidd * (2.5 + y);  // sion
    sion[1] =-KBOL * dxnidd / rho;        // dsion dd
    sion[2] = 1.5 * KBOL * dxnidd / temp; // dsion dt
    sion[3] =-y * KBOL * dxnidd * ytot;   // dsion da
    sion[4] = 0.0;                        // dsion dz
  } //helm_eos_ion

  /////////////////////////////////////////////////////////////////////////////
  // ELECTRON-POSITRON SECTION
  static void
  helm_eos_ele(const double rho, const double temp,
      double pele[5], double eele[5], double sele[5],
      double etaele[5], double xne[5],
      const struct helm_eos_cache & cache,
      const bool extrapolate = true) {
    const double ye   = cache.ye;                            // electron number fraction
    const double ytot = cache.ytot;
    const double din  = cache.din;
    int iat, jat;                                             // temperature and density indices in the table
    double fi[36];                                            // cache for the table values
    double dth, dt2, dti, dt2i, dd, dd2, ddi;                 // temperature and density deltas
    double xt, xd, mxt, mxd;                                  // various differences
    double si0t, si1t, si2t, si0mt, si1mt, si2mt;             // the six temperature basis functions
    double si0d, si1d, si2d, si0md, si1md, si2md;             // the six density basis functions
    double dsi0t, dsi1t, dsi2t, dsi0mt, dsi1mt, dsi2mt;       // derivatives of the weight functions, wrt temp
    double dsi0d, dsi1d, dsi2d, dsi0md, dsi1md, dsi2md;       // derivatives of the weight functions, wrt dens
    double ddsi0t, ddsi1t, ddsi2t, ddsi0mt, ddsi1mt, ddsi2mt; // second derivatives
    double free_en, df_d, df_t, df_tt, df_dt;                 // free energy and its derivatives
    double x, s;                                              // scratch variables
    double l10temp = log10(temp);

    jat = (l10temp - tab_ltemp_min)/tab_ltemp_delta; // hash locate temperature and density
    iat = (cache.ldin - tab_lrho_min)/tab_lrho_delta;
    if (extrapolate && (iat < tab_extrapolation_margin_irho ||
                        jat < tab_extrapolation_margin_itemp)) {
//printf("{iat, jat} = {% 3d, % 3d}\n", iat, jat);
        helm_eos_ele_offtab(rho, temp, pele, eele, sele, etaele, xne, cache);
        return;
    }

    //if (rho < tab_rho_min*(1 - HELM_EOS_EPS) || rho > tab_rho_max*(1 + HELM_EOS_EPS)) {
    if (rho > tab_rho_max*(1 + HELM_EOS_EPS)) {
      log_one(error) << "density (" << rho << ") out of table "
                     << "[" << tab_rho_min << ":" << tab_rho_max << "]" << std::endl;
      delete helm_eos_table_ptr;
      MPI_Abort(MPI_COMM_WORLD, -1);
    }
    //if (temp < tab_temp_min*(1 - HELM_EOS_EPS) || temp > tab_temp_max*(1 + HELM_EOS_EPS)) {
    if (temp > tab_temp_max*(1 + HELM_EOS_EPS)) {
      log_one(error) << "temperature (" << temp << ") out of table ["
                     << tab_temp_min << ":" << tab_temp_max << "]" << std::endl;
      delete helm_eos_table_ptr;
      MPI_Abort(MPI_COMM_WORLD, -1);
    }


    if (jat < 0 || jat >= tab_ntemp) {
      log_one(error) << "temperature (" << temp << ") off table ["
                     << 0 << ":" << tab_ntemp << "]" << std::endl;
      delete helm_eos_table_ptr;
      MPI_Abort(MPI_COMM_WORLD, -1);
    }
    if (iat < 0 || iat >= tab_nrho) {
      log_one(error) << "density (" << rho << ") off table ["
                     << 0 << ":" << tab_nrho << "]" << std::endl;
      delete helm_eos_table_ptr;
      MPI_Abort(MPI_COMM_WORLD, -1);
    }

    // compute temperature and density deltas
    dth  = helm_eos_table_ptr->temp[jat+1] - helm_eos_table_ptr->temp[jat];
    dt2  = dth * dth;
    dti  = 1.0 / dth;
    dt2i = 1.0 / dt2;
    dd   = helm_eos_table_ptr->rho[iat+1] - helm_eos_table_ptr->rho[iat];
    dd2  = dd * dd;
    ddi  = 1.0 / dd;

    // access the table locations only once
    fi[0]  = helm_eos_table_ptr->f[iat][jat];
    fi[1]  = helm_eos_table_ptr->f[iat+1][jat];
    fi[2]  = helm_eos_table_ptr->f[iat][jat+1];
    fi[3]  = helm_eos_table_ptr->f[iat+1][jat+1];
    fi[4]  = helm_eos_table_ptr->ft[iat][jat];
    fi[5]  = helm_eos_table_ptr->ft[iat+1][jat];
    fi[6]  = helm_eos_table_ptr->ft[iat][jat+1];
    fi[7]  = helm_eos_table_ptr->ft[iat+1][jat+1];
    fi[8]  = helm_eos_table_ptr->ftt[iat][jat];
    fi[9]  = helm_eos_table_ptr->ftt[iat+1][jat];
    fi[10] = helm_eos_table_ptr->ftt[iat][jat+1];
    fi[11] = helm_eos_table_ptr->ftt[iat+1][jat+1];
    fi[12] = helm_eos_table_ptr->fd[iat][jat];
    fi[13] = helm_eos_table_ptr->fd[iat+1][jat];
    fi[14] = helm_eos_table_ptr->fd[iat][jat+1];
    fi[15] = helm_eos_table_ptr->fd[iat+1][jat+1];
    fi[16] = helm_eos_table_ptr->fdd[iat][jat];
    fi[17] = helm_eos_table_ptr->fdd[iat+1][jat];
    fi[18] = helm_eos_table_ptr->fdd[iat][jat+1];
    fi[19] = helm_eos_table_ptr->fdd[iat+1][jat+1];
    fi[20] = helm_eos_table_ptr->fdt[iat][jat];
    fi[21] = helm_eos_table_ptr->fdt[iat+1][jat];
    fi[22] = helm_eos_table_ptr->fdt[iat][jat+1];
    fi[23] = helm_eos_table_ptr->fdt[iat+1][jat+1];
    fi[24] = helm_eos_table_ptr->fddt[iat][jat];
    fi[25] = helm_eos_table_ptr->fddt[iat+1][jat];
    fi[26] = helm_eos_table_ptr->fddt[iat][jat+1];
    fi[27] = helm_eos_table_ptr->fddt[iat+1][jat+1];
    fi[28] = helm_eos_table_ptr->fdtt[iat][jat];
    fi[29] = helm_eos_table_ptr->fdtt[iat+1][jat];
    fi[30] = helm_eos_table_ptr->fdtt[iat][jat+1];
    fi[31] = helm_eos_table_ptr->fdtt[iat+1][jat+1];
    fi[32] = helm_eos_table_ptr->fddtt[iat][jat];
    fi[33] = helm_eos_table_ptr->fddtt[iat+1][jat];
    fi[34] = helm_eos_table_ptr->fddtt[iat][jat+1];
    fi[35] = helm_eos_table_ptr->fddtt[iat+1][jat+1];

    // various differences
    xt  = fmax(0.0, (temp - helm_eos_table_ptr->temp[jat]) * dti);
    xd  = fmax(0.0, (din - helm_eos_table_ptr->rho[iat]) * ddi);
    mxt = 1.0 - xt;
    mxd = 1.0 - xd;

    // the six density and six temperature basis functions
    si0t  = psi0(xt);
    si1t  = psi1(xt)*dth;
    si2t  = psi2(xt)*dt2;
    si0mt = psi0(mxt);
    si1mt =-psi1(mxt)*dth;
    si2mt = psi2(mxt)*dt2;
    si0d  = psi0(xd);
    si1d  = psi1(xd)*dd;
    si2d  = psi2(xd)*dd2;
    si0md = psi0(mxd);
    si1md =-psi1(mxd)*dd;
    si2md = psi2(mxd)*dd2;

    // derivatives of the weight functions
    dsi0t  = dpsi0(xt)*dti;
    dsi1t  = dpsi1(xt);
    dsi2t  = dpsi2(xt)*dth;
    dsi0mt =-dpsi0(mxt)*dti;
    dsi1mt = dpsi1(mxt);
    dsi2mt =-dpsi2(mxt)*dth;
    dsi0d  = dpsi0(xd)*ddi;
    dsi1d  = dpsi1(xd);
    dsi2d  = dpsi2(xd)*dd;
    dsi0md =-dpsi0(mxd)*ddi;
    dsi1md = dpsi1(mxd);
    dsi2md =-dpsi2(mxd)*dd;

    // second derivatives of the weight functions
    ddsi0t  = ddpsi0(xt)*dt2i;
    ddsi1t  = ddpsi1(xt)*dti;
    ddsi2t  = ddpsi2(xt);
    ddsi0mt = ddpsi0(mxt)*dt2i;
    ddsi1mt =-ddpsi1(mxt)*dti;
    ddsi2mt = ddpsi2(mxt);

    // the free energy
    free_en = h5xy(fi,   si0t,   si1t,   si2t,   si0mt,   si1mt,   si2mt,   si0d,   si1d,   si2d,   si0md,   si1md,   si2md);
    df_t    = h5xy(fi,  dsi0t,  dsi1t,  dsi2t,  dsi0mt,  dsi1mt,  dsi2mt,   si0d,   si1d,   si2d,   si0md,   si1md,   si2md); // derivative with respect to temperature
    df_tt   = h5xy(fi, ddsi0t, ddsi1t, ddsi2t, ddsi0mt, ddsi1mt, ddsi2mt,   si0d,   si1d,   si2d,   si0md,   si1md,   si2md); // derivative with respect to temperature**2
    df_d    = h5xy(fi,   si0t,   si1t,   si2t,   si0mt,   si1mt,   si2mt,  dsi0d,  dsi1d,  dsi2d,  dsi0md,  dsi1md,  dsi2md); // derivative with respect to density
//    df_dd   = h5xy(fi,   si0t,   si1t,   si2t,   si0mt,   si1mt,   si2mt, ddsi0d, ddsi1d, ddsi2d, ddsi0md, ddsi1md, ddsi2md); // derivative with respect to density**2
    df_dt   = h5xy(fi,  dsi0t,  dsi1t,  dsi2t,  dsi0mt,  dsi1mt,  dsi2mt,  dsi0d,  dsi1d,  dsi2d,  dsi0md,  dsi1md,  dsi2md); // derivative with respect to temperature and density

    // now get the pressure derivative with density, chemical potential, and
    // electron positron number densities
    // get the interpolation weight functions
    si0t  = xpsi0(xt);
    si1t  = xpsi1(xt)*dth;
    si0mt = xpsi0(mxt);
    si1mt =-xpsi1(mxt)*dth;
    si0d  = xpsi0(xd);
    si1d  = xpsi1(xd)*dd;
    si0md = xpsi0(mxd);
    si1md =-xpsi1(mxd)*dd;

    // derivatives of weight functions
    dsi0t  = xdpsi0(xt)*dti;
    dsi1t  = xdpsi1(xt);
    dsi0mt =-xdpsi0(mxt)*dti;
    dsi1mt = xdpsi1(mxt);
    dsi0d  = xdpsi0(xd)*ddi;
    dsi1d  = xdpsi1(xd);
    dsi0md =-xdpsi0(mxd)*ddi;
    dsi1md = xdpsi1(mxd);

    // look in the pressure derivative only once
    fi[0]  = helm_eos_table_ptr->dpdf[iat][jat];
    fi[1]  = helm_eos_table_ptr->dpdf[iat+1][jat];
    fi[2]  = helm_eos_table_ptr->dpdf[iat][jat+1];
    fi[3]  = helm_eos_table_ptr->dpdf[iat+1][jat+1];
    fi[4]  = helm_eos_table_ptr->dpdft[iat][jat];
    fi[5]  = helm_eos_table_ptr->dpdft[iat+1][jat];
    fi[6]  = helm_eos_table_ptr->dpdft[iat][jat+1];
    fi[7]  = helm_eos_table_ptr->dpdft[iat+1][jat+1];
    fi[8]  = helm_eos_table_ptr->dpdfd[iat][jat];
    fi[9]  = helm_eos_table_ptr->dpdfd[iat+1][jat];
    fi[10] = helm_eos_table_ptr->dpdfd[iat][jat+1];
    fi[11] = helm_eos_table_ptr->dpdfd[iat+1][jat+1];
    fi[12] = helm_eos_table_ptr->dpdfdt[iat][jat];
    fi[13] = helm_eos_table_ptr->dpdfdt[iat+1][jat];
    fi[14] = helm_eos_table_ptr->dpdfdt[iat][jat+1];
    fi[15] = helm_eos_table_ptr->dpdfdt[iat+1][jat+1];

    // pressure derivative with density
    pele[1] = fmax(0.0, ye * h3xy(fi, si0t, si1t, si0mt, si1mt, si0d, si1d, si0md, si1md));

    // look in the electron chemical potential table only once
    fi[0]  = helm_eos_table_ptr->ef[iat][jat];
    fi[1]  = helm_eos_table_ptr->ef[iat+1][jat];
    fi[2]  = helm_eos_table_ptr->ef[iat][jat+1];
    fi[3]  = helm_eos_table_ptr->ef[iat+1][jat+1];
    fi[4]  = helm_eos_table_ptr->eft[iat][jat];
    fi[5]  = helm_eos_table_ptr->eft[iat+1][jat];
    fi[6]  = helm_eos_table_ptr->eft[iat][jat+1];
    fi[7]  = helm_eos_table_ptr->eft[iat+1][jat+1];
    fi[8]  = helm_eos_table_ptr->efd[iat][jat];
    fi[9]  = helm_eos_table_ptr->efd[iat+1][jat];
    fi[10] = helm_eos_table_ptr->efd[iat][jat+1];
    fi[11] = helm_eos_table_ptr->efd[iat+1][jat+1];
    fi[12] = helm_eos_table_ptr->efdt[iat][jat];
    fi[13] = helm_eos_table_ptr->efdt[iat+1][jat];
    fi[14] = helm_eos_table_ptr->efdt[iat][jat+1];
    fi[15] = helm_eos_table_ptr->efdt[iat+1][jat+1];

    // electron chemical potential etaele
    etaele[0] = h3xy(fi,  si0t,  si1t,  si0mt,  si1mt,  si0d,  si1d,  si0md,  si1md);
    x         = h3xy(fi,  si0t,  si1t,  si0mt,  si1mt, dsi0d, dsi1d, dsi0md, dsi1md); // derivative with respect to density
    etaele[1] = ye * x;
    etaele[2] = h3xy(fi, dsi0t, dsi1t, dsi0mt, dsi1mt,  si0d,  si1d,  si0md,  si1md); // derivative with respect to temperature
    etaele[3] =-x * din * ytot; // derivative with respect to abar and zbar
    etaele[4] = x * rho * ytot;

    // look in the number density table only once
    fi[0]  = helm_eos_table_ptr->xf[iat][jat];
    fi[1]  = helm_eos_table_ptr->xf[iat+1][jat];
    fi[2]  = helm_eos_table_ptr->xf[iat][jat+1];
    fi[3]  = helm_eos_table_ptr->xf[iat+1][jat+1];
    fi[4]  = helm_eos_table_ptr->xft[iat][jat];
    fi[5]  = helm_eos_table_ptr->xft[iat+1][jat];
    fi[6]  = helm_eos_table_ptr->xft[iat][jat+1];
    fi[7]  = helm_eos_table_ptr->xft[iat+1][jat+1];
    fi[8]  = helm_eos_table_ptr->xfd[iat][jat];
    fi[9]  = helm_eos_table_ptr->xfd[iat+1][jat];
    fi[10] = helm_eos_table_ptr->xfd[iat][jat+1];
    fi[11] = helm_eos_table_ptr->xfd[iat+1][jat+1];
    fi[12] = helm_eos_table_ptr->xfdt[iat][jat];
    fi[13] = helm_eos_table_ptr->xfdt[iat+1][jat];
    fi[14] = helm_eos_table_ptr->xfdt[iat][jat+1];
    fi[15] = helm_eos_table_ptr->xfdt[iat+1][jat+1];

    // electron + positron number densities
    xne[0] =          h3xy(fi,  si0t,  si1t,  si0mt,  si1mt,  si0d,  si1d,  si0md,  si1md);
    x      = fmax(0.0,h3xy(fi,  si0t,  si1t,  si0mt,  si1mt, dsi0d, dsi1d, dsi0md, dsi1md)); // derivative with respect to density
    xne[1] = ye * x;
    xne[2] =          h3xy(fi, dsi0t, dsi1t, dsi0mt, dsi1mt,  si0d,  si1d,  si0md,  si1md); // derivative with respect to temperature
    xne[3] =-x * din * ytot; // derivative with respect to abar and zbar
    xne[4] = x * rho * ytot;

    // Below are the desired electron-positron thermodynamic quantities

    // dpepdd at high temperatures and low densities is below the
    //   floating point limit of the subtraction of two large terms.
    // since dpresdd doesn't enter the maxwell relations at all, use the
    //   bicubic interpolation done above instead of this one
    x       = din * din;
    pele[0] = x * df_d;
    pele[2] = x * df_dt;
    //pele[1]  = ye * (x * df_dd + 2.0 * din * df_d);
    s       = pele[1] / ye - 2.0 * din * df_d;
    pele[3] = -ytot * (2.0 * pele[0] + s * din);
    pele[4] = rho * ytot * (2.0 * din * df_d  +  s);

    x       = ye * ye;
    sele[0] =-df_t * ye;
    sele[2] =-df_tt * ye;
    sele[1] =-df_dt * x;
    sele[3] = ytot * (ye * df_dt * din - sele[0]);
    sele[4] =-ytot * (ye * df_dt * rho  + df_t);

    eele[0] = ye * free_en + temp * sele[0];
    eele[2] = temp * sele[2];
    eele[1] = x * df_d + temp * sele[1];
    eele[3] =-ye * ytot * (free_en +  df_d * din) + temp * sele[3];
    eele[4] = ytot * (free_en + ye * df_d * rho) + temp * sele[4];

  } //helm_eos_ele

  /////////////////////////////////////////////////////////////////////////////
  // ELECTRON-POSITRON SECTION
  static void
  helm_eos_ele_offtab(const double rho, const double temp,
      double pele[5], double eele[5], double sele[5],
      double etaele[5], double xne[5],
      const struct helm_eos_cache & cache) {
    const double KBOL = phys::kB
                      / phys::cfactor_from_<param::cgs_units>::energy
                      * phys::cfactor_from_<param::cgs_units>::temperature;
    const double AVO = phys::NAvo;
    const double ye   = cache.ye;                            // electron number fraction
    const double ytot = cache.ytot;
    const double din  = cache.din;
    const double dxnidd = cache.dxnidd;
    const double ldin = cache.ldin;
    const double ltemp = log10(temp);

    double edeg = (ldin < eint_ele_deg_thr1) ? eint_ele_deg_coef*rho
                                             : exp10(eint_ele_deg(ldin));
    double epos = (ltemp < eint_ele_deg_thr2) ? 0.
                                              : exp10(eint_ele_ep(ltemp)
                                                - (ldin + 12.));
    double edeg_dd = (ldin < eint_ele_deg_thr1) ? eint_ele_deg_coef
                                                : dedd_ele_deg(ldin)*ye;
    double epos_dt = (ltemp < eint_ele_deg_thr2) ? 0. : dedt_ele_ep(ltemp);

    eele[0] = ye*(edeg + epos + 1.5*AVO*KBOL*temp);
    eele[1] = ye*(edeg_dd - epos/rho);
    eele[2] = epos_dt*1e-12/rho + ye*1.5*AVO*KBOL;

    double dpdd_ig = AVO*KBOL*temp*ye;
    double dpdt_ig = AVO*KBOL*din;
    double pres_ig = dpdd_ig*rho;
    double pres_deg= (ldin < -5.) ? pres_ele_deg_coef*exp10(ldin*5./3)
                                  : exp10(pres_ele_deg(ldin));
    double pres_ep = (ltemp < eint_ele_deg_thr2) ? 0.
                                                 : exp10(pres_ele_ep(ltemp));
    double dpdd_deg= (ldin < -5.) ? 5./3.*pres_deg/rho
                                  : dpdd_ele_deg(ldin)*ye;
    double dpdt_ep = (ltemp < eint_ele_deg_thr2) ? 0. : dpdt_ele_ep(ltemp);
    pele[0] = pres_ig + pres_deg + pres_ep;
    pele[1] = dpdd_ig + dpdd_deg;
    pele[2] = dpdt_ig + dpdt_ep;
  } //helm_eos_ele_offtab

  static double
  tab_eint_min_at_rho(const double rho, const double ye) {
    const double RGAS = phys::kB
                      / phys::cfactor_from_<param::cgs_units>::energy
                      * phys::cfactor_from_<param::cgs_units>::temperature
                      * phys::NAvo;

    double ldin = log10(rho*ye);
    double edeg = (ldin < eint_ele_deg_thr1) ? eint_ele_deg_coef*rho
                                             : exp10(eint_ele_deg(ldin));

    return ye*(edeg + 1.5*RGAS*temp_extrapolation_margin);
  } // tab_eint_min_at_rho

  /////////////////////////////////////////////////////////////////////////////
  // BUTTERWORTH LOW-PASS AND FIRST DERIVATIVE BY SAM JONES : "beware the butterbomb"
  struct Filter{
    double g;    // gain
    double dgdf; // d(gain)/d(frequency)
  };

  /////////////////////////////////////////////////////////////////////////////
  // BUTTERWORTH FILTER INPUT: FREQUENCY, CENTRAL FREQUENCY AND ORDER OF
  // FILTER (freq, cfreq, n), RETURNS: GAIN AND d(GAIN)/d(FREQUENCY)
  static void
  butterworth(const double freq, const double cfreq, const int n,
      struct Filter & result) {
    result.g    = 1.0/(1 + gsl_pow_int((freq / cfreq),2 * n));
    result.dgdf = -gsl_pow_2(result.g)
                * 2*n/gsl_pow_int(cfreq,2*n)*gsl_pow_int(freq,(2*n - 1));
  }

  /////////////////////////////////////////////////////////////////////////////
  // COULOMB CORRECTIONS SECTION
  static void
  helm_eos_cou(const double rho, const double temp,
      double pcoul[5], double ecoul[5], double scoul[5],
      const struct helm_eos_cache & cache) {

#   ifdef HELM_EOS_NO_COULOMB_SCREENING
    pcoul[0] = pcoul[1] = pcoul[2] = pcoul[3] = pcoul[4] = 0.;
    ecoul[0] = ecoul[1] = ecoul[2] = ecoul[3] = ecoul[4] = 0.;
    scoul[0] = scoul[1] = scoul[2] = scoul[3] = scoul[4] = 0.;
    return;
#   endif

    // fitting parameters
    const double
      a1 = -0.898004,
      b1 = 0.96786,
      c1 = 0.220703,
      d1 = -0.86097,
      e1 = 2.5269,
      a2 = 0.29561,
      b2 = 1.9885,
      c2 = 0.288675;
    //const double ye = zbar / abar; // electron number fraction
    const double KBOL = phys::kB
                      / phys::cfactor_from_<param::cgs_units>::energy
                      * phys::cfactor_from_<param::cgs_units>::temperature;
    const double HPL  = phys::hplanck
                      / phys::cfactor_from_<param::cgs_units>::energy
                      / phys::cfactor_from_<param::cgs_units>::time;
    const double AVO = phys::NAvo;
    double EE = phys::qe / phys::cfactor_from_<param::cgs_units>::charge;
    const double ytot = cache.ytot;
    const double kt = KBOL * temp;
    const double abar = cache.abar, zbar = cache.zbar;
    double xni = cache.xni, dxnidd = cache.dxnidd, dxnida = cache.dxnida;
    double pion, dpiondd, dpiondt, dpionda, dpiondz;
    double s, dsdd, dsda;
    double lami, lamidd, lamida; // average ion serperation
    double plasg, plasgdd, plasgda, plasgdt, plasgdz; // plasma coupling parameter
    struct Filter tfilter, dfilter;
    double gain, dgaindt, dgaindd;

    pion    = xni * kt;
    dpiondd = dxnidd * kt;
    dpiondt = xni * KBOL;
    dpionda = dxnida * kt;
    dpiondz = 0.0;

    s    = 4.0 / 3.0 * M_PI * xni;
    dsdd = 4.0 / 3.0 * M_PI * dxnidd;
    dsda = 4.0 / 3.0 * M_PI * dxnida;

    lami   = 1.0 / cbrt(s);
    lamidd =-lami / 3.0 * dsdd / s;
    lamida =-lami / 3.0 * dsda / s;

    plasg   = SQ(EE * zbar) / (kt * lami);
    plasgdd =-plasg / lami * lamidd;
    plasgda =-plasg / lami * lamida;
    plasgdt =-plasg / temp;
    plasgdz = 2.0 * plasg / zbar;

    if (plasg >= 1.0) {
      double x = sqrt(sqrt(plasg));
      double y = AVO * KBOL * ytot;
      ecoul[0] = y * temp * (a1 * plasg + b1 * x + c1 / x + d1);
      pcoul[0] = rho * ecoul[0]  / 3.0;
      scoul[0] =-y * (3.0 * b1 * x - 5.0 * c1 / x + d1 * (log(plasg) - 1.0) - e1);

      y        = AVO * KBOL * temp * ytot * (a1 + 0.25 / plasg * (b1 * x - c1 / x));
      ecoul[1] = y * plasgdd;
      ecoul[2] = y * plasgdt + ecoul[0] / temp;
      ecoul[3] = y * plasgda - ecoul[0] / abar;
      ecoul[4] = y * plasgdz;

      y        = rho / 3.0;
      pcoul[1] = ecoul[0] / 3.0 + y * ecoul[1];
      pcoul[2] = y * ecoul[2];
      pcoul[3] = y * ecoul[3];
      pcoul[4] = y * ecoul[4];

      y        =-AVO * KBOL / (abar / plasg) * (0.75 * b1 * x + 1.25 * c1 / x + d1);
      scoul[1] = y * plasgdd;
      scoul[2] = y * plasgdt;
      scoul[3] = y * plasgda - scoul[0] / abar;
      scoul[4] = y * plasgdz;
    }
    else {
      const double x = plasg * sqrt(plasg);
      const double y = pow(plasg, b2);
      const double z = c2 * x - a2 / 3.0 * y;

      pcoul[0] =-pion * z;
      ecoul[0] = 3.0 * pcoul[0] / rho;
      scoul[0] =-AVO * KBOL / abar * (c2 * x - a2 * (b2 - 1.0) / b2 * y);

      s        = 1.5 * c2 * x / plasg - a2 * b2 / 3.0 * y / plasg;
      pcoul[1] =-dpiondd * z - pion * s * plasgdd;
      pcoul[2] =-dpiondt * z - pion * s * plasgdt;
      pcoul[3] =-dpionda * z - pion * s * plasgda;
      pcoul[4] =-dpiondz * z - pion * s * plasgdz;

      s        = 3.0 / rho;
      ecoul[1] = s * pcoul[1] - ecoul[0] / rho;
      ecoul[2] = s * pcoul[2];
      ecoul[3] = s * pcoul[3];
      ecoul[4] = s * pcoul[4];

      s        =-AVO * KBOL / (abar * plasg) * (1.5 * c2 * x - a2 * (b2 - 1.0) * y);
      scoul[1] = s * plasgdd;
      scoul[2] = s * plasgdt;
      scoul[3] = s * plasgda - scoul[0] / abar;
      scoul[4] = s * plasgdz;
    }

    // butterworth bomb proofing by Sam Jones : "beware the butterbomb"
    butterworth(log10(temp) - 4.5, 3.0, 12, tfilter);
    butterworth(log10(rho)  + 1.0, 6.0, 12, dfilter);

    // derivatives (and conversion from logarithmic derivative)
    gain    = (1.0 - tfilter.g * dfilter.g);
    dgaindt =-dfilter.g * tfilter.dgdf / temp / M_LN10;
    dgaindd =-tfilter.g * dfilter.dgdf / rho / M_LN10;

    // derivatives via chain rule
    pcoul[1] = gain * pcoul[1] + pcoul[0] * dgaindd;
    pcoul[2] = gain * pcoul[2] + pcoul[0] * dgaindt;
    pcoul[3] = gain * pcoul[3];
    pcoul[4] = gain * pcoul[4];

    ecoul[1] = gain * ecoul[1] + ecoul[0] * dgaindd;
    ecoul[2] = gain * ecoul[2] + ecoul[0] * dgaindt;
    ecoul[3] = gain * ecoul[3];
    ecoul[4] = gain * ecoul[4];

    scoul[1] = gain * scoul[1] + scoul[0] * dgaindd;
    scoul[2] = gain * scoul[2] + scoul[0] * dgaindt;
    scoul[3] = gain * scoul[3];
    scoul[4] = gain * scoul[4];

    // straight up gain
    pcoul[0] = pcoul[0] * gain;
    ecoul[0] = ecoul[0] * gain;
    scoul[0] = scoul[0] * gain;

  } //helm_eos_cou

  static void
  eos_helm_init(const char* datafile) {
    FILE *file;

    file = fopen(datafile, "r");
    std::ifstream infile(datafile);
    if (!infile.good()) {
      log_one(error) << "error opening EOS table file" << std::endl << "the filname was " << datafile << std::endl;
      MPI_Abort(MPI_COMM_WORLD, -1);
    }

    //helm_eos_table_ptr = safe_malloc<helm_eos_table>(sizeof(struct helm_eos_table));
    helm_eos_table_ptr = new helm_eos_table();

    // READ THE HELMHOLTZ FREE ENERGY TABLE AND ITS DERIVATIVES
    for (int j = 0; j < tab_ntemp; j++) {
      helm_eos_table_ptr->temp[j] = exp10(tab_ltemp_min + j*tab_ltemp_delta);
      for (int i = 0; i < tab_nrho; i++) {
        helm_eos_table_ptr->rho[i] = exp10(tab_lrho_min + i*tab_lrho_delta);
        if (fscanf(file, "%lf %lf %lf %lf %lf %lf %lf %lf %lf",
          &helm_eos_table_ptr->f[i][j], &helm_eos_table_ptr->fd[i][j], &helm_eos_table_ptr->ft[i][j],
          &helm_eos_table_ptr->fdd[i][j], &helm_eos_table_ptr->ftt[i][j], &helm_eos_table_ptr->fdt[i][j],
          &helm_eos_table_ptr->fddt[i][j], &helm_eos_table_ptr->fdtt[i][j], &helm_eos_table_ptr->fddtt[i][j]) != 9) {
  	      log_one(error) << "error reading the Helmholtz free energy table at i = " << i << ", j = " << j << std::endl;
  	      fclose(file);
  	      delete (helm_eos_table_ptr);
          MPI_Abort(MPI_COMM_WORLD, -1);
        }
      }
    }

    // READ THE PRESSURE TABLE AND ITS DERIVATIVES
    for (int j = 0; j < tab_ntemp; j++) {
      for (int i = 0; i < tab_nrho; i++) {
        if (fscanf(file, "%lf %lf %lf %lf",
          &helm_eos_table_ptr->dpdf[i][j], &helm_eos_table_ptr->dpdfd[i][j],
          &helm_eos_table_ptr->dpdft[i][j], &helm_eos_table_ptr->dpdfdt[i][j]) != 4) {
          log_one(error) << "error reading the pressure derivative table at i = " << i << ", j = " << j << std::endl;
  	      fclose(file);
  	      delete helm_eos_table_ptr;
          MPI_Abort(MPI_COMM_WORLD, -1);
        }
      }
    }

    // READ THE ELECTRON CHEMICAL POTENTIAL TABLE AND ITS DERIVATIVES
    for (int j = 0; j < tab_ntemp; j++) {
      for (int i = 0; i < tab_nrho; i++) {
        if (fscanf(file, "%lf %lf %lf %lf",
          &helm_eos_table_ptr->ef[i][j], &helm_eos_table_ptr->efd[i][j],
          &helm_eos_table_ptr->eft[i][j], &helm_eos_table_ptr->efdt[i][j]) != 4) {
          log_one(error) << "error reading the electron chemical potential table at i = " << i << ", j = " << j << std::endl;
  	      fclose(file);
  	      delete helm_eos_table_ptr;
          MPI_Abort(MPI_COMM_WORLD, -1);
        }
      }
    }

    // READ THE NUMBER DENSITY TABLE AND ITS DERIVATIVES
    for (int j = 0; j < tab_ntemp; j++) {
      for (int i = 0; i < tab_nrho; i++) {
        if (fscanf(file, "%lf %lf %lf %lf",
          &helm_eos_table_ptr->xf[i][j], &helm_eos_table_ptr->xfd[i][j],
          &helm_eos_table_ptr->xft[i][j], &helm_eos_table_ptr->xfdt[i][j]) != 4) {
  	      log_one(error) << "error reading the number density table at i = " << i << ", j = " << j << std::endl;
  	      fclose(file);
  	      delete helm_eos_table_ptr;
          MPI_Abort(MPI_COMM_WORLD, -1);
        }
      }
    }

    fclose(file);

    // Initialize extrapolation function
    eos_helm_extrapolation_init();

  } //eos_helm_init

  static void
  eos_helm_extrapolation_init() {
    // Initialize boundary interpolating functions
    const double KBOL = phys::kB
                      / phys::cfactor_from_<param::cgs_units>::energy
                      * phys::cfactor_from_<param::cgs_units>::temperature;
    const double AVO = phys::NAvo;

    struct helm_eos_cache cache;
    double pele[5]={0}, eele[5]={0}, sele[5]={0}, etaele[5]={0}, xne[5]={0};
    double rho_margin = exp10(tab_lrho_min
                      + tab_extrapolation_margin_irho*tab_lrho_delta);
    double temp_margin = exp10(tab_ltemp_min
                       + tab_extrapolation_margin_itemp*tab_ltemp_delta);
    // log10(rho_margin)  = -11.5
    // log10(temp_margin) =   3.5

    // low temperature boundary
    double eint_ltemp[tab_nrho], pres_ltemp[tab_nrho];

    for (int i = 0; i < tab_nrho; i++) {
      double rho = exp10(tab_lrho_min + i*tab_lrho_delta);
      helm_eos_update_cache(rho, 1., 1., cache);
      double lrho = cache.ldin;
      helm_eos_ele(rho,temp_margin,pele,eele,sele,etaele,xne,cache,false);
      eint_ltemp[i] = eele[0];

      double p_lin = pres_ele_deg_coef2*exp10(lrho);
      pres_ltemp[i] = log10(std::abs(pele[0] - p_lin));
    } // i

    double eint0 = eint_ltemp[10];
    for (int i = 0; i < tab_nrho; i++) {
      eint_ltemp[i] = log10(std::abs(eint_ltemp[i] - eint0));
    } // i

    // density derivative
    double dedd_ltemp[tab_nrho], dpdd_ltemp[tab_nrho];
    for (int i = 0; i < tab_nrho; i++) {
      int i1 = std::max(0, i - 1);
      int i2 = std::min(tab_nrho - 1, i + 1);
      double lrho = tab_lrho_min + i*tab_lrho_delta;
      dedd_ltemp[i] = (eint_ltemp[i2] - eint_ltemp[i1])
                    / ((i2 - i1)*tab_lrho_delta)
                    * exp10(eint_ltemp[i] - lrho);
      dpdd_ltemp[i] = (pres_ltemp[i2] - pres_ltemp[i1])
                    / ((i2 - i1)*tab_lrho_delta)
                    * exp10(pres_ltemp[i] - lrho);
    } // i

    eint_ele_deg.set_data(tab_lrho_min, tab_lrho_max, eint_ltemp, tab_nrho);
    dedd_ele_deg.set_data(tab_lrho_min, tab_lrho_max, dedd_ltemp, tab_nrho);
    pres_ele_deg.set_data(tab_lrho_min, tab_lrho_max, pres_ltemp, tab_nrho);
    dpdd_ele_deg.set_data(tab_lrho_min, tab_lrho_max, dpdd_ltemp, tab_nrho);

    // low temperature boundary
    double eint_lrho[tab_ntemp], pres_lrho[tab_ntemp];

    helm_eos_update_cache(tab_rho_min, 1., 1., cache);
    for (int i = 0; i < tab_ntemp; i++) {
      double temp = exp10(tab_ltemp_min + i*tab_ltemp_delta);
      helm_eos_ele(rho_margin,temp,pele,eele,sele,etaele,xne,cache,false);
      eint_lrho[i] = log10(std::abs(eele[0] - eint0 * temp / 1e3));
      double pgas = AVO*KBOL*temp*rho_margin;
      pres_lrho[i] = log10(std::abs(pele[0] - pgas));
    } // i

    // temperature derivative
    double dedt_lrho[tab_ntemp], dpdt_lrho[tab_nrho];
    for (int i = 0; i < tab_ntemp; i++) {
      int i1 = std::max(0, i - 1);
      int i2 = std::min(tab_ntemp - 1, i + 1);
      double ltemp = tab_ltemp_min + i*tab_ltemp_delta;
      dedt_lrho[i] = (eint_lrho[i2] - eint_lrho[i1])
                   / ((i2 - i1)*tab_ltemp_delta)
                   * exp10(eint_lrho[i] - ltemp);
      dpdt_lrho[i] = (pres_lrho[i2] - pres_lrho[i1])
                   / ((i2 - i1)*tab_ltemp_delta)
                   * exp10(pres_lrho[i] - ltemp);
    } // i
    eint_ele_ep.set_data(tab_ltemp_min, tab_ltemp_max, eint_lrho, tab_ntemp);
    dedt_ele_ep.set_data(tab_ltemp_min, tab_ltemp_max, dedt_lrho, tab_ntemp);
    pres_ele_ep.set_data(tab_ltemp_min, tab_ltemp_max, pres_lrho, tab_ntemp);
    dpdt_ele_ep.set_data(tab_ltemp_min, tab_ltemp_max, dpdt_lrho, tab_ntemp);
    //for (int i=0; i<tab_ntemp; ++i) {
    //  double ltemp = tab_ltemp_min + i*tab_ltemp_delta;
    //  printf ("%15.12f  %24.15e\n", ltemp, eint_ele_ep(ltemp));
    //}

  } //eos_helm_extrapolation_init

  /////////////////////////////////////////////////////////////////////////////
  // FREE TABLE MEMORY
  static void eos_deinit() {
    if (helm_eos_table_ptr == NULL) return;
    delete helm_eos_table_ptr;
  }

  /////////////////////////////////////////////////////////////////////////////
  // UPDATE CACHE FOR OFTEN USED QUANTITIES
  static void
  helm_eos_update_cache(double rho, double abar, double zbar,
      struct helm_eos_cache & cache) {
    const double AVO = phys::NAvo;
    cache.abar   = abar;
    cache.zbar   = zbar;
    cache.ytot   = 1.0 / abar;
    cache.ye     = zbar * cache.ytot;
    cache.din    = rho * cache.ye;
    cache.ldin   = log10(cache.din);
    cache.xni    = AVO * cache.ytot * rho;
    cache.dxnidd = AVO * cache.ytot;
    cache.dxnida =-cache.xni * cache.ytot;
  } // helm_eos_update_cache

  /////////////////////////////////////////////////////////////////////////////
  // GETTING PRESSURE, SOUNDSPEED, AND E_INT FROM RHO AND TEMP
  static void
  helm_eos_given_rho_t(body & b) {
    double rho  = b.getDensity();
    double temp = b.getTemperature();
    double m    = b.mass();
    double abar = b.getAbar();
    double zbar = abar*b.getElectronfraction();
    struct helm_eos_cache cache;

    helm_eos_update_cache(rho, abar, zbar, cache);
    double prad[5] = {0}, pion[5] = {0}, pele[5] = {0}, pcou[5] = {0};
    double erad[5] = {0}, eion[5] = {0}, eele[5] = {0}, ecou[5] = {0};
    double srad[5] = {0}, sion[5] = {0}, sele[5] = {0}, scou[5] = {0};
    double etaele[5] = {0}, xne[5] = {0};
    double p[5] = {0}, e[5] = {0}, s[5] = {0};

    helm_eos_rad(rho, temp, prad, erad, srad);
    helm_eos_ion(rho, temp, pion, eion, sion, cache);
    helm_eos_ele(rho, temp, pele, eele, sele, etaele, xne, cache);
    helm_eos_cou(rho, temp, pcou, ecou, scou, cache);
    for (int i = 0; i < 5; i++) {
      p[i] = prad[i] + pion[i] + pele[i] + pcou[i];
      e[i] = erad[i] + eion[i] + eele[i] + ecou[i];
      s[i] = srad[i] + sion[i] + sele[i] + scou[i];
    }

    double cv       = e[DTEMP];
    double chit     = temp/p[VALUE]*p[DTEMP];
    double chid     = p[DRHO]*rho/p[VALUE];
    double x        = p[VALUE]/rho * chit/(temp*cv);
    double gamma_1  = chit*x + chid;
    const double C_LIGHT_CGS
                 = phys::clight
                 / phys::cfactor_from_<param::cgs_units>::velocity;
    double sound = C_LIGHT_CGS*sqrt(gamma_1
                   /(1 + (e[VALUE] + SQ(C_LIGHT_CGS))*rho/p[VALUE]));

    b.setPressure(p[VALUE]);
    b.setSoundspeed(sound);
    b.setInternalenergy(e[VALUE]);
    b.setEntropy(s[VALUE]);
  } //helm_eos_given_rho_t

  /////////////////////////////////////////////////////////////////////////////
  // GETTING PRESSURE, INT_E, SOUNDSPEED, AND TEMPERATURE FROM RHO AND S
  static void
  helm_eos_given_rho_s(body & b) {
    double ent = b.getEntropy(); // entropy used for convergence
    double m = b.mass(),         rho = b.getDensity(),
        abar = b.getAbar(),     zbar = abar*b.getElectronfraction(),
       _temp = b.getTemperature(); // temperature first guess
    int iter;                   // number of Newton-Raphson iterations
    double _s = 0.0, _dt = 0.0;

    struct helm_eos_cache cache;
    double e[5] = {0}, p[5] = {0}, s[5] = {0};
    double temp;

    helm_eos_update_cache(rho, abar, zbar, cache);
    double prad[5] = {0}, pion[5] = {0}, pele[5] = {0}, pcou[5] = {0};
    double erad[5] = {0}, eion[5] = {0}, eele[5] = {0}, ecou[5] = {0};
    double srad[5] = {0}, sion[5] = {0}, sele[5] = {0}, scou[5] = {0};
    double etaele[5] = {0}, xne[5] = {0};
    if (param::convergence_method == param::newton_raphson) {
      for (iter = 0; iter < HELM_EOS_MAXITER; iter++) {
        helm_eos_rad(rho, _temp, prad, erad, srad);
        helm_eos_ion(rho, _temp, pion, eion, sion, cache);
        helm_eos_ele(rho, _temp, pele, eele, sele, etaele, xne, cache);
        helm_eos_cou(rho, _temp, pcou, ecou, scou, cache);
        temp = _temp;
        for (int i = 0; i < 5; i++) {
          p[i] = prad[i] + pion[i] + pele[i] + pcou[i];
          e[i] = erad[i] + eion[i] + eele[i] + ecou[i];
          s[i] = srad[i] + sion[i] + sele[i] + scou[i];
        }
        _s = s[VALUE];
        if ( fabs(_s - ent) < fabs(HELM_EOS_EPS*ent) ) break;
        _dt = -(_s - ent) / s[DTEMP];
        if ((_temp + _dt) <= tab_temp_min) {
          _temp = tab_temp_min;
          break;
        }
        if ((_temp + _dt) >= tab_temp_max) {
          _temp = tab_temp_max;
          break;
        }
        _temp += _dt;
      }
      if (iter >= HELM_EOS_MAXITER) {
        log_one(error) << "Newton-Raphson in function did not converge." << std::endl;
        delete helm_eos_table_ptr;
        MPI_Abort(MPI_COMM_WORLD, -1);
      }
    } else if (param::convergence_method == param::bisection) {
      double T_a = tab_temp_min, T_b = tab_temp_max;
      double T_c = 0.0;
      double s_c = 0.0;
      helm_eos_rad(rho, T_a, prad, erad, srad);
      helm_eos_ion(rho, T_a, pion, eion, sion, cache);
      helm_eos_ele(rho, T_a, pele, eele, sele, etaele, xne, cache);
      helm_eos_cou(rho, T_a, pcou, ecou, scou, cache);
      double s_a = srad[0] + sion[0] + sele[0] + scou[0];
      helm_eos_rad(rho, T_b, prad, erad, srad);
      helm_eos_ion(rho, T_b, pion, eion, sion, cache);
      helm_eos_ele(rho, T_b, pele, eele, sele, etaele, xne, cache);
      helm_eos_cou(rho, T_b, pcou, ecou, scou, cache);
      double s_b = srad[0] + sion[0] + sele[0] + scou[0];
      for (iter = 0; iter < HELM_EOS_MAXITER; iter++) {
        T_c = 0.5*(T_a + T_b);
        helm_eos_rad(rho, T_c, prad, erad, srad);
        helm_eos_ion(rho, T_c, pion, eion, sion, cache);
        helm_eos_ele(rho, T_c, pele, eele, sele, etaele, xne, cache);
        helm_eos_cou(rho, T_c, pcou, ecou, scou, cache);
        s_c = srad[0] + sion[0] + sele[0] + scou[0];
        if (  (fabs(s_c - ent) < fabs(HELM_EOS_EPS*ent))
           or (fabs(T_b - T_a) < HELM_EOS_EPS) ) {
          temp = T_c;
          for (int i = 0; i < 5; i++) {
            p[i] = prad[i] + pion[i] + pele[i] + pcou[i];
            e[i] = erad[i] + eion[i] + eele[i] + ecou[i];
            s[i] = srad[i] + sion[i] + sele[i] + scou[i];
          }
          break;
        }
        if (std::copysign(1.0, s_c - ent) == std::copysign(1.0,s_a - ent)) {
          T_a = T_c;
          s_a = s_c;
        } else {
          T_b = T_c;
          s_b = s_c;
        }
      }
      if (iter >= HELM_EOS_MAXITER) {
        log_one(error) << "Bisection in function did not converge." << std::endl;
        delete helm_eos_table_ptr;
        MPI_Abort(MPI_COMM_WORLD, -1);
      }
    }
    double cv       = e[DTEMP];
    double chit     = _temp/p[VALUE]*p[DTEMP];
    double chid     = p[DRHO]*rho/p[VALUE];
    double x        = p[VALUE]/rho * chit/(temp*cv);
    double gamma_1  = chit*x + chid;
    const double C_LIGHT_CGS
                 = phys::clight
                 / phys::cfactor_from_<param::cgs_units>::velocity;
    double sound = C_LIGHT_CGS*sqrt(gamma_1
                   / (1 + (e[VALUE] + SQ(C_LIGHT_CGS))*rho/p[VALUE]));

    b.setPressure(p[VALUE]);
    b.setInternalenergy(e[VALUE]);
    b.setSoundspeed(sound);
    b.setTemperature(temp);
    // b.setGamma(cp/cv);
  } //helm_eos_given_rho_s

  /////////////////////////////////////////////////////////////////////////////
  // ENTROPY: RADIATION SECTION
  static double
  entropy_helm_eos_rad(const double rho, const double temp) {
    const double SIG_OVER_C = phys::sB * phys::clight
                            / phys::cfactor_from_<param::cgs_units>::energy
                            * phys::cfactor_from_<param::cgs_units>::volume;
    return 16.*SIG_OVER_C*CU(temp)/(3.*rho);
  } // entropy_helm_eos_rad

  /////////////////////////////////////////////////////////////////////////////
  // ENTROPY: ION SECTION
  static double
  entropy_helm_eos_ion(const double rho, const double temp,
      const struct helm_eos_cache & cache) {
    const double abar = cache.abar;
    const double ytot = cache.ytot;
    const double xni  = cache.xni;
    const double KBOL = phys::kB
                      / phys::cfactor_from_<param::cgs_units>::energy
                      * phys::cfactor_from_<param::cgs_units>::temperature;
    const double AMU  = phys::amu
                      / phys::cfactor_from_<param::cgs_units>::mass;
    const double HPL  = phys::hplanck
                      / phys::cfactor_from_<param::cgs_units>::energy
                      / phys::cfactor_from_<param::cgs_units>::time;
    const double AVO = phys::NAvo;
    //const double y = cache.ywot + cache.lswot15 + 1.5 * ltemp;
    const double s = (2.0 * M_PI * AMU * KBOL) / SQ(HPL) * temp;
    const double y = log(abar * abar * sqrt(abar) / (rho * AVO) * s * sqrt(s));

    return KBOL*(2.5*xni/rho + AVO*ytot*y);        // sion
  } //helm_eos_ion

  /////////////////////////////////////////////////////////////////////////////
  // ENTROPY: ELECTRON-POSITRON SECTION
  static double
  entropy_helm_eos_ele(const double rho, const double temp,
      const struct helm_eos_cache & cache) {
    const double ye   = cache.ye;                        // electron number fraction
    const double ytot = cache.ytot;
    const double din  = cache.din;
    int iat, jat;                                        // temperature and density indices in the table
    double fi[36];                                       // cache for the table values
    double dth, dti, dd, dd2, ddi;                       // temperature and density deltas
    double xt, xd, mxt, mxd;                             // various differences
    double si0d, si1d, si2d, si0md, si1md, si2md;        // the six density basis functions
    double dsi0t, dsi1t, dsi2t, dsi0mt, dsi1mt, dsi2mt;  // derivatives of the weight functions, wrt temp
    double df_t;                                         // free energy and its derivatives
    double l10temp = log(temp) * (1.0 / log(10.0));

    if (rho < tab_rho_min || rho > tab_rho_max) {
      log_one(error) << "density (" << rho << ") out of table "
                     << "[" << tab_rho_min << ":" << tab_rho_max << "]" << std::endl;
      delete helm_eos_table_ptr;
      MPI_Abort(MPI_COMM_WORLD, -1);
    }
    if (temp < tab_temp_min || temp > tab_temp_max) {
      log_one(error) << "temperature (" << temp << ") out of table ["
                     << tab_temp_min << ":" << tab_temp_max << "]" << std::endl;
      delete helm_eos_table_ptr;
      MPI_Abort(MPI_COMM_WORLD, -1);
    }

    jat = (l10temp - tab_ltemp_min)/tab_ltemp_delta; // hash locate temperature and density
    iat = (cache.ldin - tab_lrho_min)/tab_lrho_delta;

    if (jat < 0 || jat >= tab_ntemp) {
      log_one(error) << "temperature (" << temp << ") off table ["
                     << 0 << ":" << tab_ntemp << "]" << std::endl;
      delete helm_eos_table_ptr;
      MPI_Abort(MPI_COMM_WORLD, -1);
    }
    if (iat < 0 || iat >= tab_nrho) {
      log_one(error) << "density (" << rho << ") produced index ("
                     << iat << ") off table [" << 0 << ":" << tab_nrho << "]"
                     << std::endl;
      delete helm_eos_table_ptr;
      MPI_Abort(MPI_COMM_WORLD, -1);
    }

    // compute temperature and density deltas
    dth  = helm_eos_table_ptr->temp[jat+1] - helm_eos_table_ptr->temp[jat];
    dti  = 1.0 / dth;
    dd   = helm_eos_table_ptr->rho[iat+1] - helm_eos_table_ptr->rho[iat];
    ddi  = 1.0 / dd;
    dd2  = dd*dd;

    // access the table locations only once
    fi[0]  = helm_eos_table_ptr->f[iat][jat];
    fi[1]  = helm_eos_table_ptr->f[iat+1][jat];
    fi[2]  = helm_eos_table_ptr->f[iat][jat+1];
    fi[3]  = helm_eos_table_ptr->f[iat+1][jat+1];
    fi[4]  = helm_eos_table_ptr->ft[iat][jat];
    fi[5]  = helm_eos_table_ptr->ft[iat+1][jat];
    fi[6]  = helm_eos_table_ptr->ft[iat][jat+1];
    fi[7]  = helm_eos_table_ptr->ft[iat+1][jat+1];
    fi[8]  = helm_eos_table_ptr->ftt[iat][jat];
    fi[9]  = helm_eos_table_ptr->ftt[iat+1][jat];

    fi[10] = helm_eos_table_ptr->ftt[iat][jat+1];
    fi[11] = helm_eos_table_ptr->ftt[iat+1][jat+1];
    fi[12] = helm_eos_table_ptr->fd[iat][jat];
    fi[13] = helm_eos_table_ptr->fd[iat+1][jat];
    fi[14] = helm_eos_table_ptr->fd[iat][jat+1];
    fi[15] = helm_eos_table_ptr->fd[iat+1][jat+1];
    fi[16] = helm_eos_table_ptr->fdd[iat][jat];
    fi[17] = helm_eos_table_ptr->fdd[iat+1][jat];
    fi[18] = helm_eos_table_ptr->fdd[iat][jat+1];
    fi[19] = helm_eos_table_ptr->fdd[iat+1][jat+1];
    fi[20] = helm_eos_table_ptr->fdt[iat][jat];
    fi[21] = helm_eos_table_ptr->fdt[iat+1][jat];
    fi[22] = helm_eos_table_ptr->fdt[iat][jat+1];
    fi[23] = helm_eos_table_ptr->fdt[iat+1][jat+1];
    fi[24] = helm_eos_table_ptr->fddt[iat][jat];
    fi[25] = helm_eos_table_ptr->fddt[iat+1][jat];
    fi[26] = helm_eos_table_ptr->fddt[iat][jat+1];
    fi[27] = helm_eos_table_ptr->fddt[iat+1][jat+1];
    fi[28] = helm_eos_table_ptr->fdtt[iat][jat];
    fi[29] = helm_eos_table_ptr->fdtt[iat+1][jat];
    fi[30] = helm_eos_table_ptr->fdtt[iat][jat+1];
    fi[31] = helm_eos_table_ptr->fdtt[iat+1][jat+1];
    fi[32] = helm_eos_table_ptr->fddtt[iat][jat];
    fi[33] = helm_eos_table_ptr->fddtt[iat+1][jat];
    fi[34] = helm_eos_table_ptr->fddtt[iat][jat+1];
    fi[35] = helm_eos_table_ptr->fddtt[iat+1][jat+1];

    // various differences
    xt  = fmax(0.0, (temp - helm_eos_table_ptr->temp[jat]) * dti);
    xd  = fmax(0.0, (din - helm_eos_table_ptr->rho[iat]) * ddi);
    mxt = 1.0 - xt;
    mxd = 1.0 - xd;

    // the six density basis functions
    si0d  = psi0(xd);
    si1d  = psi1(xd)*dd;
    si2d  = psi2(xd)*dd2;
    si0md = psi0(mxd);
    si1md =-psi1(mxd)*dd;
    si2md = psi2(mxd)*dd2;

    // derivatives of the weight functions
    dsi0t  = dpsi0(xt)*dti;
    dsi1t  = dpsi1(xt);
    dsi2t  = dpsi2(xt)*dth;
    dsi0mt =-dpsi0(mxt)*dti;
    dsi1mt = dpsi1(mxt);
    dsi2mt =-dpsi2(mxt)*dth;

    // the free energy:
    // derivative with respect to temperature
    df_t    = h5xy(fi,  dsi0t,  dsi1t,  dsi2t,  dsi0mt, dsi1mt, dsi2mt,
                         si0d,   si1d,   si2d,   si0md,  si1md,  si2md);

    return -df_t * ye;
  } // entropy_helm_eos_ele

  /////////////////////////////////////////////////////////////////////////////
  // COULOMB CORRECTIONS SECTION
  static double
  entropy_helm_eos_cou(const double rho, const double temp,
      const struct helm_eos_cache & cache) {
    // return value:
    double ent = 0.0;

    // fitting parameters
    const double
      a1 =-0.898004,
      b1 = 0.96786,
      c1 = 0.220703,
      d1 =-0.86097,
      e1 = 2.5269,
      a2 = 0.29561,
      b2 = 1.9885,
      c2 = 0.288675;

    const double KBOL = phys::kB
                      / phys::cfactor_from_<param::cgs_units>::energy
                      * phys::cfactor_from_<param::cgs_units>::temperature;
    const double AVO = phys::NAvo;
    double ytot  = cache.ytot;
    double kt    = KBOL * temp;
    double abar  = cache.abar, zbar = cache.zbar;
    double xni   = cache.xni;
    double s     = 4.0 / 3.0 * M_PI * xni;
    double lami  = 1.0 / cbrt(s);
    double EE = phys::qe / phys::cfactor_from_<param::cgs_units>::charge;
    double plasg = SQ(EE * zbar) / (kt * lami);

    if (plasg >= 1.0) {
      double x = sqrt(sqrt(plasg));
      double y = AVO * KBOL * ytot;
      ent =-y * (3.0 * b1 * x - 5.0 * c1 / x + d1 * (log(plasg) - 1.0) - e1);
    }
    else {
      double x = plasg * sqrt(plasg);
      double y = pow(plasg, b2);
      double z = c2 * x - a2 / 3.0 * y;

      ent =-AVO * KBOL / abar * (c2 * x - a2 * (b2 - 1.0) / b2 * y);
    }

    // butterworth bomb proofing by Sam Jones : "beware the butterbomb"
    struct Filter tfilter, dfilter;
    butterworth(log10(temp) - 4.5, 3.0, 12, tfilter);
    butterworth(log10(rho)  + 1.0, 6.0, 12, dfilter);

    // derivatives (and conversion from logarithmic derivative)
    double gain = (1.0 - tfilter.g * dfilter.g);
    return ent * gain;

  } // entropy_helm_eos_cou

  /**
  * @brief      Helper function: computes internal energy and derivatives
  *             given density and temperature
  *
  * @param      rho:    density
  * @param      temp:   temperature
  * @param      eint:   output array for int. energy and derivatives
  */
  static void
  get_eint_given_rho_temp (const double rho, const double temp, double eint[5],
      const struct helm_eos_cache & cache) {

    double prad[5] = {0}, pion[5] = {0}, pele[5] = {0}, pcou[5] = {0};
    double erad[5] = {0}, eion[5] = {0}, eele[5] = {0}, ecou[5] = {0};
    double srad[5] = {0}, sion[5] = {0}, sele[5] = {0}, scou[5] = {0};
    double etaele[5] = {0}, xne[5] = {0};

    helm_eos_rad(rho, temp, prad, erad, srad);
    helm_eos_ion(rho, temp, pion, eion, sion, cache);
    helm_eos_ele(rho, temp, pele, eele, sele, etaele, xne, cache);
    helm_eos_cou(rho, temp, pcou, ecou, scou, cache);

    for (int i = 0; i < 5; ++i)
      eint[i] = erad[i] + eion[i] + eele[i] + ecou[i];

  } // get_eint_given_rho_temp


}; // class eos_t<param::eos_helmholtz>

#if eos_type == eos_helmholtz
  get_quantity_t eos_t<param::eos_helmholtz>::get_dpdrho_at_temp = nullptr;
#endif

eos_t<param::eos_helmholtz>::helm_eos_table* eos_t<param::eos_helmholtz>::helm_eos_table_ptr = nullptr;
struct interpolating_function_1d eos_t<param::eos_helmholtz>::eint_ele_deg{};
struct interpolating_function_1d eos_t<param::eos_helmholtz>::dedd_ele_deg{};
struct interpolating_function_1d eos_t<param::eos_helmholtz>::eint_ele_ep{};
struct interpolating_function_1d eos_t<param::eos_helmholtz>::dedt_ele_ep{};

struct interpolating_function_1d eos_t<param::eos_helmholtz>::pres_ele_deg{};
struct interpolating_function_1d eos_t<param::eos_helmholtz>::dpdd_ele_deg{};
struct interpolating_function_1d eos_t<param::eos_helmholtz>::pres_ele_ep{};
struct interpolating_function_1d eos_t<param::eos_helmholtz>::dpdt_ele_ep{};

} // namespace eos
