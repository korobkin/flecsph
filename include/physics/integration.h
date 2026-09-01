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
    double rp = flecsph::magnitude(pos);
    double vr = vel[0];
    double t1 = physics::totaltime_prev + influx::extraction_radius
              / vr*(influx::extraction_radius/rp - 1.);
    double rn = influx::extraction_radius
              / (1. - vr*(physics::totaltime - t1)/influx::extraction_radius);
    pos *= rn/rp;
    source.set_coordinates(pos);
    if (rn > influx::extraction_radius) {
      // Set up particles for hydro evolution once past the extraction sphere
      //printf("for particle %08d t1=%12.5f, totaltime=%12.5f, and totaltime_prev=%12.5f\n",source.id(),t1,physics::totaltime,physics::totaltime_prev);
      source.set_state(NONE);
      const double psi = (rn - influx::extraction_radius)/(rn-rp);
      double dt = psi*physics::dt;
      double dt1 = physics::totaltime + physics::dt - t1;
      if (rn < influx::extraction_radius + param::flow_velocity*physics::dt*C_LIGHT_CGS) {
        pos *= influx::extraction_radius/rn;
      }
      const double
        R2 = pos[0]*pos[0] + pos[1]*pos[1],
        R = sqrt(R2 + 1e-12),
        r2 = R2 + pos[2]*pos[2],
        r = sqrt(r2 + 1e-12),
        cos_phi = pos[0]/R,    sin_phi = pos[1]/R,
        cos_tht = pos[2]/r,    sin_tht = R/r,
        vth = vel[1],          vphi = vel[2];
      double phi = atan2(sin_phi, cos_phi);
      if(phi < 0){phi += 2*M_PI;}
      if (rn < influx::extraction_radius + param::flow_velocity*physics::dt*C_LIGHT_CGS) {
        influx::grid_data_point_t gp = influx::linear_interpolator(physics::totaltime_prev + dt, atan2(sin_tht,cos_tht), phi);
        printf("for particle %08d dt=%12.5f, dt1=%12.5f, totaltime=%12.5f, and totaltime_prev=%12.5f\n",source.id(),dt,dt1,physics::totaltime,physics::totaltime_prev);
        vr = gp.vr*C_LIGHT_CGS;
      }
      //printf("for particle %08d vr_old=%12.5e, vr_new=%12.5e\n",source.id(),vr/C_LIGHT_CGS, gp.vr*C_LIGHT_CGS);
      vel[0] = vr*sin_tht*cos_phi + vth*cos_tht*cos_phi + vphi*cos_phi,
      vel[1] = vr*sin_tht*sin_phi + vth*cos_tht*sin_phi + vphi*sin_phi,
      vel[2] = vr*cos_tht - vth*sin_tht;
      source.setVelocity(vel);
      source.setVelocityhalf(vel);
      // take step with proper velocities for remainder of time step
      if (rn < influx::extraction_radius + param::flow_velocity*physics::dt*C_LIGHT_CGS) {
        pos[0] += vel[0]*dt;
        pos[1] += vel[1]*dt;
        pos[2] += vel[2]*dt;
      }
      //printf("for particle %08d rn/Rex=%12.5e, rp=%12.5e, r_n=%12.5e, r=%12.5e\n",source.id(),rn/influx::extraction_radius, rp, rn, r);
      //printf("%08d %12.5e %12.5e %12.5e\n",source.id(),sqrt(pos[0]*pos[0]+pos[1]*pos[1]+pos[2]*pos[2]),gp.rho,vr/C_LIGHT_CGS);
      double rr = sqrt(pos[0]*pos[0]+pos[1]*pos[1]+pos[2]*pos[2]);
      vr = pos[0]*vel[0]+pos[1]*vel[1]+pos[2]*vel[2]/rr;
      printf("%08d %12.5e %12.5e %12.5e\n",source.id(),rr,source.getDensity(),vr);
      //printf("for particle %08d vr_diff=%12.5e\n",source.id(), vr/C_LIGHT_CGS - gp.vr);
      // TODO:: Add a flag here to only run this when using input flux or when we want this calculated
      // Calculate difference in flux and gradient of flux at the boundary
     }
  }
  else {
    source.set_coordinates(
        source.coordinates() + physics::dt * source.getVelocity());
  }
}

}; // namespace integration

#endif // _integration_h_
