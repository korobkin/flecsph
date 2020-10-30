/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2018 Triad National Security, LLC
 * All rights reserved.
 * --------------------------------------------------------------------------~*/

/******************************************************************************
 *                                                                            *
 * EOS_UTILS.H *
 *                                                                            *
 * Global macros, utilities including root dinfer, definitiions, etc          *
 *                                                                            *
 ******************************************************************************/

#ifndef _eos_utils_h_
#define _eos_utils_h_

#include <gsl/gsl_eigen.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_vector.h>

#include "eos_consts.h"
#include "params.h"

namespace eos{
  template<param::eos_type_keyword EOS_TYPE>
  class eos_t{
  };

  // main eos function type
  typedef void (*compute_quantity_t)(body &);
  typedef double (*get_quantity_t)(const body &);

} // namespace eos

// Passive variables (if present)
#define PASSIVE_START (NVAR_BASE)
#define PASSIVE_STOP (NVAR_BASE + NVAR_PASSIVE)
#define PASSTYPE_INTRINSIC (0)
#define PASSTYPE_NUMBER (1)
#define YE (PASSIVE_START)

// EOS
#define EOS_TYPE_GAMMA (0)
#define EOS_TYPE_POLYTROPE (1)
#define EOS_TYPE_TABLE (2)
#define EOS_NUM_EXTRA (0)
#define EOS_LRHO (0)
#define EOS_LT (1)
#define EOS_YE (2)
// mass fractions
#define NUM_MASS_FRACTIONS (4)
#define MF_XA (0)
#define MF_XH (1)
#define MF_XN (2)
#define MF_XP (3)

// Fixup parameters
// may only apply for EOS GAMMA
constexpr double RHOMINLIMIT = 1.e-17;
constexpr double UUMINLIMIT = 1.e-20;
constexpr double RHOMIN = 1.e-5;
constexpr double UUMIN = 1.e-8;
constexpr double BSQORHOMAX = 50.;
constexpr double BSQOUMAX = 2500.;
constexpr double RHOEPS = 2.0;
constexpr double UORHOMAX = 50.;

// Numerical convenience to represent a small (<< 1) non-zero quantity
constexpr double SMALL = 1.e-20;

// Loop over primitive variables
#define PLOOP for(int ip = 0; ip < NVAR; ip++)
#define BASELOOP for(int ip = 0; ip < NVAR_BASE; ip++)

// Loop over extra variables
// TODO: Figure out how to make this conditionally defined.
#define EOS_ELOOP for(int e = 0; e < EOS_NUM_EXTRA; e++)

// ----------------------------------------------------------------------
// Function defs
// TODO : Make it correct order. We define above here because of ordering
double EOS_Poly_pressure_rho0_u(double rho, double u, double K, double Gam);
double EOS_Poly_pressure_rho0_w(double rho, double w, double K, double Gam);
double EOS_Poly_enthalpy_rho0_u(double rho, double u, double K, double Gam);
double EOS_Poly_entropy_rho0_u(double rho, double u, double K, double Gam);
double EOS_Poly_sound_speed_rho0_u(double rho, double u, double K, double Gam);
double EOS_Poly_rho_floor(double scale, double bsq);
double EOS_Poly_u_floor(double scale, double bsq);
void EOS_Poly_set_floors(double scale,
  double rho,
  double u,
  double bsq,
  double * rhoflr,
  double * uflr);
double EOS_Poly_adiabatic_constant(double rho, double u, double K, double Gam);
double EOS_Poly_temperature(double rho, double u, double K, double Gam);
double EOS_Poly_u_press(double press, double rho, double K, double Gam);
double EOS_Poly_Theta_unit();
void * safe_malloc(int size);
void safe_system(const char * command);
void safe_fscanf(FILE * stream, const char * format, ...);
double find_min(const double * array, int size);
double find_max(const double * array, int size);
int find_index(double value, const double * array, int size);
double interp_1d(double x,
  const double xmin,
  const double xmax,
  const int imin,
  const int imax,
  const double * tab_x,
  const double * tab_y);
void set_units();

// Structs
// ----------------------------------------------------------------------
struct of_adiabat {
  double s, ye;
  double lrho_min, lrho_max;
  int imin, imax;
  double hm1_min, hm1_max;
  double * lT;
};

struct of_tablebounds {
  int Nrho, NT, NYe;
  double lrho_min, lrho_max, dlrho;
  double lT_min, lT_max, dlT;
  double Ye_min, Ye_max, dYe;
};

// ----------------------------------------------------------------------

// Important global variables
// ----------------------------------------------------------------------
struct GV {
  static double rho_poly_thresh;
  static double poly_K, poly_gam;
  static double Reh;
  static double Risco;
  static double mbh, Mbh, L_unit, T_unit, M_unit, RHO_unit, U_unit, B_unit;
  static double TEMP_unit;
};
// Default values for Global Variables
double GV::B_unit = 1;
double GV::L_unit = 1;
double GV::M_unit = 1;
double GV::T_unit = 1;
double GV::U_unit = 1;
double GV::RHO_unit = 1;
double GV::TEMP_unit = 1;

// Public function APIs
// ----------------------------------------------------------------------

// Adding polytrope eos for fallbakcing
double
EOS_Poly_pressure_rho0_u(double rho, double u, double K, double Gam) {
  rho = fabs(rho + SMALL);
  return K * pow(rho, Gam);
}

double
EOS_Poly_pressure_rho0_w(double rho, double w, double K, double Gam) {
  rho = fabs(rho + SMALL);
  return K * pow(rho, Gam);
}
double
EOS_Poly_enthalpy_rho0_u(double rho, double u, double K, double Gam) {
  double P = EOS_Poly_pressure_rho0_u(rho, u, K, Gam);
  return rho + u + P;
}

double
EOS_Poly_entropy_rho0_u(double rho, double u, double K, double Gam) {
  // TODO: units?
  return K;
}

double
EOS_Poly_sound_speed_rho0_u(double rho, double u, double K, double Gam) {
  rho = std::max(0.0, rho);
  u = std::max(0.0, u);
  double rhogam = K * pow(rho, Gam - 1);
  double cs2_num = Gam * rhogam;
  double cs2_den = 1 + u + rhogam;
  double cs2 = cs2_num / cs2_den;
  return sqrt(cs2);
}

double
EOS_Poly_rho_floor(double scale, double bsq) {
  double rhoflr = RHOMIN * scale;
  rhoflr = std::max(rhoflr, RHOMINLIMIT);
  rhoflr = std::max(rhoflr, bsq / BSQORHOMAX);

  return rhoflr;
}

double
EOS_Poly_u_floor(double scale, double bsq) {
  double uflr = UUMIN * scale;
  uflr = std::max(uflr, UUMINLIMIT);
  uflr = std::max(uflr, bsq / BSQOUMAX);

  return uflr;
}

void
EOS_Poly_set_floors(double scale,
  double rho,
  double u,
  double bsq,
  double * rhoflr,
  double * uflr) {
  *rhoflr = EOS_Poly_rho_floor(scale, bsq);
  *uflr = EOS_Poly_u_floor(scale, bsq);
}

double
EOS_Poly_adiabatic_constant(double rho, double u, double K, double Gam) {
  return K / (Gam - 1.);
}

double
EOS_Poly_temperature(double rho, double u, double K, double Gam) {
  // TODO: is this right?
  return 0.0; // Polytrope is at zero internal energy/temperature
}

double
EOS_Poly_u_press(double press, double rho, double K, double Gam) {
  return 0.0; // pressure and internal energy independent
}

double
EOS_Poly_Theta_unit() {
  return MP / ME;
}

// Utilities
template<typename T>
T *
safe_malloc(int size) {
  size *= sizeof(T);
  // malloc(0) may or may not return NULL, depending on compiler.
  if(size == 0)
    return NULL;

  T * A = static_cast<T *>(malloc(size));
  if(A == NULL) {
    fprintf(stderr, "Failed to malloc\n");
    exit(-1);
  }
  return A;
}

// Error-handling wrappers for standard C functions
void
safe_system(const char * command) {
  int systemReturn = system(command);
  if(systemReturn == -1) {
    fprintf(stderr, "system() call %s failed! Exiting!\n", command);
    exit(-1);
  }
}

void
safe_fscanf(FILE * stream, const char * format, ...) {
  va_list args;
  va_start(args, format);
  int vfscanfReturn = vfscanf(stream, format, args);
  va_end(args);
  if(vfscanfReturn == -1) {
    fprintf(stderr, "fscanf() call failed! Exiting!\n");
    exit(-1);
  }
}

double
find_min(const double * array, int size) {
  double min = INFINITY;
  for(int i = 0; i < size; i++) {
    if(array[i] < min)
      min = array[i];
  }
  return min;
}

double
find_max(const double * array, int size) {
  double max = -INFINITY;
  for(int i = 0; i < size; i++) {
    if(array[i] > max)
      max = array[i];
  }
  return max;
}

int
find_index(double value, const double * array, int size) {
  for(int i = 0; i < size; i++) {
    if(array[i] >= value)
      return i;
  }
  return -1;
}

// A very general 1D linear interpolator
double
interp_1d(double x,
  const double xmin,
  const double xmax,
  const int imin,
  const int imax,
  const double * tab_x,
  const double * tab_y) {
  if(x < xmin)
    x = xmin;
  if(x > xmax)
    x = xmax - SMALL;
  const int nx = imax - imin;
  const double dx = (xmax - xmin) / (nx - 1);
  const int ix = imin + (x - xmin) / dx;
  const double delx = (x - tab_x[ix]) / dx;
  const double out = (1 - delx) * tab_y[ix] + delx * tab_y[ix + 1];
  /*
  // DEBUGGING
  if (isnan(out) || x > xmax || x < xmin) {
    fprintf(stderr,"[interp_1d]: out is NaN!\n");
    fprintf(stderr,"\t\tx           = %g\n",x);
    fprintf(stderr,"\t\timin        = %d\n",imin);
    fprintf(stderr,"\t\timax        = %d\n",imax);
    fprintf(stderr,"\t\txmin        = %e\n",xmin);
    fprintf(stderr,"\t\txmax        = %e\n",xmax);
    fprintf(stderr,"\t\tnx          = %d\n",nx);
    fprintf(stderr,"\t\tdx          = %e\n",dx);
    fprintf(stderr,"\t\tix          = %d\n",ix);
    fprintf(stderr,"\t\tdelx        = %e\n",delx);
    fprintf(stderr,"\t\ttab_x[ix]   = %e\n",tab_x[ix]);
    fprintf(stderr,"\t\ttab_x[ix+1] = %e\n",tab_x[ix+1]);
    fprintf(stderr,"\t\ttab_y[ix]   = %e\n",tab_y[ix]);
    fprintf(stderr,"\t\ttab_y[ix+1] = %e\n",tab_y[ix+1]);
    fprintf(stderr,"\n");
    exit(1);
  }
  */
  return out;
}

void
set_units() {
  GV::T_unit = GV::L_unit / C_LIGHT_CGS;
  GV::RHO_unit = GV::M_unit * pow(GV::L_unit, -3.);
  GV::U_unit = GV::RHO_unit * C_LIGHT_CGS * C_LIGHT_CGS;
  GV::B_unit = C_LIGHT_CGS * sqrt(4. * M_PI * GV::RHO_unit);
  GV::TEMP_unit = KBOL / MEV; // temp(MeV)/GV::TEMP_UNIT = K
}

#endif
