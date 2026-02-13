/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2017 Triad National Security, LLC
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
 * @file integration.h
 * @author Julien Loiseau
 * @date October 2018
 * @brief Integration methods
 */

#ifndef _integration_h_
#define _integration_h_

#include <vector>

#include "default_physics.h"
#include "params.h"
#include "influx.h"

namespace integration {
using namespace param;

/**
 * @brief      v -> v12
 *
 * @param      srch  The source's body holder
 */
void
save_velocityhalf(body & source) {
  source.setVelocityhalf(source.getVelocity());
}

/**
 * @brief      Leapfrog: kick velocity
 *             v^{n+1/2} = v^{n} + (dv/dt)^n * dt/2
 *             or
 *             v^{n+1} = v^{n+1/2} + (dv/dt)^n * dt/2
 *
 * @param      srch  The source's body holder
 */
void
leapfrog_kick_v(body & source) {
  if (enable_inflow and source.state() == INACTIVE) return;
  source.setVelocity(
    source.getVelocity() +
    0.5 * physics::dt * (source.getAcceleration() + source.getGAcceleration()));
}

/**
 * @brief      Leapfrog: kick internal energy
 *             u^{n+1/2} = u^{n} + (du/dt)^n * dt/2
 *             or
 *             u^{n+1} = u^{n+1/2} + (du/dt)^n * dt/2
 *
 * @param      particle
 */
void
leapfrog_kick_u(body & particle) {
  if (enable_inflow and particle.state() == INACTIVE) return;
  const double du = 0.5 * physics::dt * particle.getDudt();
  const double eint = particle.getInternalenergy();

  if (eint + du < 0.0)
    particle.setInternalenergy(eint*exp(du/eint));
  else
    particle.setInternalenergy(eint + du);
    
}

/**
 * @brief      Leapfrog: kick thermokinetic or total energy
 *             e^{n+1/2} = e^{n} + (de/dt)^n * dt/2
 *             or
 *             e^{n+1} = e^{n+1/2} + (de/dt)^n * dt/2
 *
 * @param      srch  The source's body holder
 */
void
leapfrog_kick_e(body & source) {
  if (enable_inflow and source.state() == INACTIVE) return;
  source.setTotalenergy(
    source.getTotalenergy() + 0.5 * physics::dt * source.getDedt());
}

/**
 * @brief      Leapfrog: drift
 *             r^{n+1} = r^{n} + v^{n+1/2} * dt
 *
 * @param      srch  The source's body holder
 */
void
leapfrog_drift(body & source) {
  if (enable_inflow and (source.state() == INACTIVE)) {
    point_t pos = source.coordinates();
    point_t vel = source.getVelocity();
    double rp = flecsi::magnitude(pos);
    double vr = vel[0];
    double t1 = physics::totaltime_prev + influx::extraction_radius
              / vr*(influx::extraction_radius/rp - 1.);
    double rn = influx::extraction_radius
              / (1. - vr*(physics::totaltime - t1)/influx::extraction_radius);
    pos *= rn/rp;
    source.set_coordinates(pos);
    if (rn > influx::extraction_radius) {
      // Set up particles for hydro evolution once past the extraction sphere
      source.set_state(NONE);
      const double 
        R2 = pos[0]*pos[0] + pos[1]*pos[1],
        R = sqrt(R2 + 1e-12),
        r2 = R2 + pos[2]*pos[2],
        r = sqrt(r2 + 1e-12),
        cos_phi = pos[0]/R,    sin_phi = pos[1]/R,
        cos_tht = pos[2]/r,    sin_tht = R/r,
        vth = vel[1],          vphi = vel[2];
      vel[0] = vr*sin_tht*cos_phi + vth*cos_tht*cos_phi + vphi*cos_phi,
      vel[1] = vr*sin_tht*sin_phi + vth*cos_tht*sin_phi + vphi*sin_phi,
      vel[2] = vr*cos_tht - vth*sin_tht;
      source.setVelocity(vel);
      source.setVelocityhalf(vel);
      // TODO:: Add a flag here to only run this when using input flux or when we want this calculated
      // Calculate difference in flux and gradient of flux at the boundary
      double phi = atan2(sin_phi, cos_phi);
      if(phi < 0){phi += 2*M_PI;}
      
      if(influx::input_flux_unfolded_in_hdf5){
         const double vel_corr = (param::flow_velocity*C_LIGHT_CGS)/vr;
         const double r_corr = (param::flow_velocity*C_LIGHT_CGS*t1)/param::sphere_radius;
         influx::grid_data_point_t gp = influx::linear_interpolator(t1, atan2(sin_tht,cos_tht), phi)*vel_corr*r_corr*r_corr;
         //influx::grid_data_point_t gp = influx::linear_interpolator(t1, atan2(sin_tht,cos_tht), phi);
         double diff = gp.rho - source.getDensity();
         printf("%08d %12.5f %12.5f %12.5e %12.5f %12.5f %12.5f\n", source.id(), physics::totaltime, vel_corr, r_corr, gp.rho, source.getDensity(), diff);
      }
      else {
         influx::grid_data_point_t gp = influx::linear_interpolator(physics::totaltime, atan2(sin_tht,cos_tht), phi);
         double diff = gp.rho - source.getDensity();
         printf("%08d %12.5f %12.5f %12.5f %12.5f\n", source.id(), physics::totaltime, gp.rho, source.getDensity(), diff);
      }
  // const point_t gradrho_prof = density_profiles::grad_rho_ndim_from_data_grid(pos);
    }
  }
  else {
    source.set_coordinates(
        source.coordinates() + physics::dt * source.getVelocity());
  }
}

}; // namespace integration

#endif // _integration_h_
