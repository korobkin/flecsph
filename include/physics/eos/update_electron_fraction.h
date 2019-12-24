#ifndef _update_electron_fraction_h_
#define _update_electron_fraction_h_

#include <iostream>
#include <hdf5.h>
#include "eos_preamble.h"
#include "eos_utils.h"

namespace update_ye {

  void read_ye_table(const char *name)
  {

  std::ifstream infile(name);  // for the following code to work, this must be a .hdf5 file
  if(!infile.good()){
    clog_one(error)<<"File "<<name<<" not found."<<std::endl;
    MPI_Finalize();
    exit(-1);
  }


  hsize_t file_grid_dims[3], file_start[3], file_count[3];
  hsize_t mem_grid_dims[3], mem_start[3];

  hid_t file_id = H5Fopen(name, H5F_ACC_RDONLY, H5P_DEFAULT); // opens existing hdf5 file
  
  hid_t dset_id = H5Dopen(file_id, "/pointsrho", H5P_DEFAULT);
  int status = H5Dread(dset_id, H5T_STD_I32LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &Nrho_ye);
  H5Dclose(dset_id);

  dset_id = H5Dopen(file_id, "/pointsT", H5P_DEFAULT);
  status = H5Dread(dset_id, H5T_STD_I32LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &NT_ye);
  H5Dclose(dset_id);

  dset_id = H5Dopen(file_id, "/pointsye_initial", H5P_DEFAULT);
  status = H5Dread(dset_id, H5T_STD_I32LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &NYe_ye);
  H5Dclose(dset_id);

  int tab_size = Nrho_ye*NT_ye*NYe_ye;

  tab_rho = safe_malloc<double>(Nrho_ye);
  tab_T = safe_malloc<double>(NT_ye);
  tab_Yeye = safe_malloc<double>(NYe_ye);

  tab_dYedt = safe_malloc<double>(tab_size);
  tab_deweakdt = safe_malloc<double>(tab_size);

  dset_id = H5Dopen(file_id, "/rho", H5P_DEFAULT);
  status = H5Dread(dset_id, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, tab_rho);
  H5Dclose(dset_id);

  dset_id = H5Dopen(file_id, "/T", H5P_DEFAULT);
  status = H5Dread(dset_id, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, tab_T);
  H5Dclose(dset_id);

  dset_id = H5Dopen(file_id, "/ye_initial", H5P_DEFAULT);
  status = H5Dread(dset_id, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, tab_Yeye);
  H5Dclose(dset_id);

  dset_id = H5Dopen(file_id, "/dyedt", H5P_DEFAULT);
  status = H5Dread(dset_id, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, tab_dYedt);
  H5Dclose(dset_id);

  dset_id = H5Dopen(file_id, "/de_weakdt", H5P_DEFAULT);
  status = H5Dread(dset_id, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, tab_deweakdt);
  H5Dclose(dset_id);

  H5Fclose(file_id);

  tab_rhoye_min = find_min(tab_rho, Nrho_ye);
  tab_rhoye_max = find_max(tab_rho, Nrho_ye);
  tab_Tye_min = find_min(tab_T, NT_ye);
  tab_Tye_max = find_max(tab_T, NT_ye);
  tab_Yeye_min = find_min(tab_Yeye, NYe_ye);
  tab_Yeye_max = find_max(tab_Yeye, NYe_ye);
  tab_dYedt_min = find_min(tab_dYedt, tab_size);
  tab_deweakdt_min = find_min(tab_deweakdt, tab_size);

  tab_drhoye = (tab_rhoye_max - tab_rhoye_min)/(Nrho_ye-1);
  tab_dTye = (tab_Tye_max - tab_Tye_min)/(NT_ye-1);
  tab_dYeye = (tab_Yeye_max - tab_Yeye_min)/(NYe_ye-1);
  }

  static double catch_rhoye(const double rho)
  {
    if (rho < tab_rhoye_min) {
      return tab_rhoye_min;
    } else if (rho > tab_rhoye_max) {
      return tab_rhoye_max - TABLE_TOL/100.;
    } else {
      return rho;
    }
  }

  static double catch_Tye(const double T)
  {
    if (T < tab_Tye_min) {
      return tab_Tye_min;
    } else if (T >= tab_Tye_max) {
      return tab_Tye_max - TABLE_TOL/100.;
    } else {
      return T;
    }
  }

  static double catch_yeye(const double ye)
  {
    if ( ye < tab_Yeye_min ) {
      return tab_Yeye_min;
    } else if ( ye >= tab_Yeye_max ) {
      return tab_Yeye_max - TABLE_TOL/100.;
    } else {
      return ye;
    }
  }

  // Interpolation
  // ----------------------------------------------------------------------
  static double YE_interp(const double rho, const double T, const double Ye,
            double* tab)
  {
    // We don't want our indices to fall off the table,
    // but we want to extrapolate appropriately
    double T2    = catch_Tye(T);
    double Ye2   = catch_yeye(Ye);
    double rho2  = catch_rhoye(rho);


    // indices
    const int irho = (rho2 - tab_rhoye_min)/tab_drhoye;
    const int iT   = (T2   - tab_Tye_min)/tab_dTye;
    const int iY   = (Ye2   - tab_Yeye_min)/tab_dYeye;

    // assumes evenly spaced table
    const double delrho = (rho - tab_rho[irho])/tab_drhoye;
    const double delT   = (T   - tab_T[iT])/tab_dTye;
    const double delY   = (Ye   - tab_Yeye[iY])/tab_dYeye;

    // one-sided, trilinear interpolation
    return (
      (1-delrho)    * (1-delT)  * (1-delY)  * tab[YE_ELEM(irho,   iT,   iY)]
      + delrho      * (1-delT)  * (1-delY)  * tab[YE_ELEM(irho+1, iT,   iY)]
      + (1-delrho)  * delT      * (1-delY)  * tab[YE_ELEM(irho,   iT+1, iY)]
      + delrho      * delT      * (1-delY)  * tab[YE_ELEM(irho+1, iT+1, iY)]
      + (1-delrho)  * (1-delT)  * delY      * tab[YE_ELEM(irho,   iT,   iY+1)]
      + delrho      * (1-delT)  * delY      * tab[YE_ELEM(irho+1, iT,   iY+1)]
      + (1-delrho)  * delT      * delY      * tab[YE_ELEM(irho,   iT+1, iY+1)]
      + delrho      * delT      * delY      * tab[YE_ELEM(irho+1, iT+1, iY+1)]
      );
  }

  // ----------------------------------------------------------------------

  double compute_deweakdt(double rho, double T, double ye_initial)
  {
    return YE_interp(rho, T, ye_initial, tab_deweakdt);
  }

  // Testing to see if code is written correctly for navigating dYe/dt table
  double get_dYedt(body& b) {
    double rho = b.getDensity();
    double T = b.getTemperature();
    double ye_init = b.getElectronfraction();
    double dYedt_find = YE_interp(log10(rho), T/1.0e9, ye_init, tab_dYedt);
    return dYedt_find;
  }
}


#endif
