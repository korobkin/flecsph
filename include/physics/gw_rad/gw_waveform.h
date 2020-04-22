/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2020 Triad National Security, LLC
 * All rights reserved.
 *~--------------------------------------------------------------------------~*/

 /*~--------------------------------------------------------------------------~*
 *
 * /@@@@@@@@  @@           @@@@@@   @@@@@@@@ @@@@@@@  @@      @@
 * /@@/////  /@@          @@////@@ @@////// /@@////@@/@@     /@@
 * /@@       /@@  @@@@@  @@    // /@@       /@@   /@@/@@     /@@
 * /@@@@@@@  /@@ @@///@@/@@       /@@@@@@@@@/@@@@@@@ /@@@@@@@@@@
 * /@@////   /@@/@@@@@@@/@@       ////////@@/@@////  /@@//////@@
 * /@@       /@@/@@//// //@@    @@       /@@/@@      /@@     /@@
 * /@@       @@@//@@@@@@ //@@@@@@  @@@@@@@@ /@@      /@@     /@@
 * //       ///  //////   //////  ////////  //       //      //
 *
 *~--------------------------------------------------------------------------~*/

/**
 * @file gw_waveform.h
 * @authore Hyun Lim
 * @date Apr 2020
 * @brief Gravitatioanl waveform extraction for Newtonian source
 *        by calculating quadrupole moments in slow motion approximation
 */

#ifndef _GW_WAVEFORM_H_
#define _GW_WAVEFORM_H_

#include "utils_gw.h"

/*
 * Extract GW information via angle-averaged strain values
 * for + and x polarization
 */

void
extract_gw_waveform(std::vector<body> &bodies) {

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
   
   //Define angle averaged value of strain:
   // <rh_+> and <rh_x>.

   double strain_hp = 0.0; // Plus polarization
   double strain_hc = 0.0; // Cross polarization

   // Ref : Zhuge et al. PRD.50.6247, 1994
   //       Blanchet. LRR-2014-2
   //       van den Broek et al. MNRAS 425, L24-L27, 2012

      // HL : Here, I propose some alternative way to compute rhx and rh+
      //      This is based on just using particles' position, velocity,
      //      and acceleration to compute quadrupole moments and its derivatives

      double qxx_dtdt = 0.0, qyy_dtdt = 0.0, qzz_dtdt = 0.0;
      double qxy_dtdt = 0.0, qxz_dtdt = 0.0, qyz_dtdt = 0.0;

      for (auto b:bodies){
        // Compute second time derivatives of each components of quadrupole moments
        qxx_dtdt = 2.0/3.0 * b.mass()
                 * (2.0*b.coordinates()[0]*b.getAcceleration()[0]
                    - b.coordinates()[1]*b.getAcceleration()[1]
                    - b.coordinates()[2]*b.getAcceleration()[2]
                    + 2.0*b.getVelocity()[0]*b.getVelocity()[0]
                    - b.getVelocity()[1]*b.getVelocity()[1]
                    - b.getVelocity()[2]*b.getVelocity()[2]
                   );

        qxy_dtdt = b.mass()
                 * (b.coordinates()[0]*b.getAcceleration()[1]
                    + b.coordinates()[1]*b.getAcceleration()[0]
                    + 2.0*b.getVelocity()[0]*b.getVelocity()[1]
                   );

        qxz_dtdt = b.mass()
                 * (b.coordinates()[0]*b.getAcceleration()[2]
                    + b.coordinates()[2]*b.getAcceleration()[0]
                    + 2.0*b.getVelocity()[0]*b.getVelocity()[2]
                   );

        qyy_dtdt = 2.0/3.0 * b.mass()
                 * (2.0*b.coordinates()[1]*b.getAcceleration()[1]
                    - b.coordinates()[0]*b.getAcceleration()[0]
                    - b.coordinates()[2]*b.getAcceleration()[2]
                    + 2.0*b.getVelocity()[1]*b.getVelocity()[1]
                    - b.getVelocity()[0]*b.getVelocity()[0]
                    - b.getVelocity()[2]*b.getVelocity()[2]
                   );

        qyz_dtdt = b.mass()
                 * (b.coordinates()[1]*b.getAcceleration()[2]
                    + b.coordinates()[2]*b.getAcceleration()[1]
                    + 2.0*b.getVelocity()[1]*b.getVelocity()[2]
                   );

        qzz_dtdt = 2.0/3.0 * b.mass()
                 * (2.0*b.coordinates()[2]*b.getAcceleration()[2]
                    - b.coordinates()[1]*b.getAcceleration()[1]
                    - b.coordinates()[0]*b.getAcceleration()[0]
                    + 2.0*b.getVelocity()[2]*b.getVelocity()[2]
                    - b.getVelocity()[1]*b.getVelocity()[1]
                    - b.getVelocity()[0]*b.getVelocity()[0]
                   );

        // Using symmetric property defining remaining values
        double qyx_dtdt = qxy_dtdt;
        double qzx_dtdt = qxz_dtdt;
        double qzy_dtdt = qyz_dtdt;
      }

     mpi_utils::reduce_sum(qxx_dtdt);
     mpi_utils::reduce_sum(qxy_dtdt);
     mpi_utils::reduce_sum(qxz_dtdt);
     mpi_utils::reduce_sum(qyy_dtdt);
     mpi_utils::reduce_sum(qyz_dtdt);
     mpi_utils::reduce_sum(qzz_dtdt);

     // Prefactors : depends on unit system

     double cl = C_LIGHT_CGS;
     double Gc = param::gravitational_constant;

     double prefac_L = Gc/(cl*cl*cl*cl*cl);
     double prefac_h = cl*cl*cl*cl/Gc;
     double prefac_h_sq = prefac_h*prefac_h;


     // Compute angle averaged strain values
     double strain_hp_sq = 0.0, strain_hc_sq = 0.0;

     strain_hp_sq = 4.0/15.0*((qxx_dtdt - qzz_dtdt)*(qxx_dtdt - qzz_dtdt)
                              + (qyy_dtdt - qzz_dtdt)*(qyy_dtdt - qzz_dtdt)
                              + qxz_dtdt*qxz_dtdt + qyz_dtdt*qyz_dtdt)
                  + 1.0/10.0*(qxx_dtdt - qyy_dtdt)*(qxx_dtdt - qyy_dtdt)
                  + 14.0/15.0*qxy_dtdt*qxy_dtdt/prefac_h_sq;

     strain_hc_sq = 1.0/6.0*(qxx_dtdt - qyy_dtdt)*(qxx_dtdt - qyy_dtdt)
                  + 2.0/3.0*qxy_dtdt*qxy_dtdt
                  + 4.0/3.0*(qxz_dtdt*qxz_dtdt + qyz_dtdt*qyz_dtdt)/prefac_h_sq;

     strain_hp = std::sqrt(strain_hp_sq);
     strain_hc = std::sqrt(strain_hc_sq);

     #if 0
     // HL : I just added this for future ref
     if(do_2p5PN){
     //TODO : add higher order PN
     }
     #endif

     //Compute gravitational luminosity

     double Lgw = 0.0;
     double ang_I3t = 0.0; //HL : will finish this formulation

     Lgw = prefac_L*ang_I3t/5.0;

     //TODO : Checking now. Define var for now to check outputting option

     // Getting output
     //TODO : Make it separate function? Move this to analysis.h?
     static bool first_time = true;
     if (param::out_scalar_every <=0 ||
         physics::iteration % param::out_scalar_every !=0)
         return;

     // output only from rank 0
     if (rank !=0) return;
     const char *filename = "gw_waveform_info.dat";

     if (first_time) {
       // Generate output header
       /* HL :  Here, I assume that we only accept 3 dimensional case.
        *       It system of dimension is less than 3, code will be stopped
        * TODO : will be generalized once we have lower dimensional case
        *         such as axisymmetry case
        */
       std::ostringstream oss_header;
       switch(gdimension){
       case 1:
         std::cerr<<"System of dimension must be 3"<<std::endl;
         assert(false);
       break;

       case 2:
         std::cerr<<"System of dimension must be 3"<<std::endl;
         assert(false);
       break;

       case 3:
       default:
         oss_header
           << "# GW Waveform Data:"<<std::endl
           << "# 1:iteration 2:time 3:timestep "
           << "# 4:rh+ 5:rhx 6:Lgw"<<std::endl;
       }

       std::ofstream out(filename);
       out << oss_header.str();
       out.close();
       first_time = false;
     }

     //TODO : Better way?
     std::ostringstream oss_data;
     oss_data << std::setw(14) << physics::iteration
        << std::setw(20) << std::scientific << std::setprecision(12)
        << physics::totaltime << std::setw(20) << physics::dt << " "
        << strain_hp << " " << strain_hc << " " << Lgw <<std::endl;

     // Open file in append mode
     std::ofstream out(filename,std::ios_base::app);
     out << oss_data.str();
     out.close();
   // If we don't include this routine, we will
   // only have zeros for strain TODO : Good?

 } //Evaluate GW waveform

#endif //_GW_WAVEFORM_H_
