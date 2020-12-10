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
 * @file apm.h
 * @date December 2020
 * @brief Interface to select various artificial pressure computation
 */

#pragma once

#include "params.h"
#include "tree.h"
#include "user.h"
#include <boost/algorithm/string.hpp>
#include <math.h>
#include <stdlib.h>
#include "lane_emden.h"
#include "density_profiles.h"
#include "kernels.h"

#define SQ(x) ((x) * (x))
#define CU(x) ((x) * (x) * (x))

namespace apm {

template<param::apm_type_keyword K>
point_t compute_apm_acc(body & particle, std::vector<body *> &nbs);

typedef point_t (*compute_apm_acc_t)(body & particle, std::vector<body *> &nbs);

/**
 * @brief  Zero APM force
 */
template<>
point_t
compute_apm_acc<param::zero_apm>(body & particle, std::vector<body *> &nbs) {

    point_t a_apm = 0.0;

    return a_apm;

}
/**
 * @brief  kilonova spherical-ejecta 
 */
template<>
point_t
compute_apm_acc<param::kn_ejecta>(body & particle, std::vector<body *> &nbs) {

    point_t a_apm = 0.0;

    return a_apm;

}

/**
 * @brief  Sharp spherical density profile
 */
template<>
point_t
compute_apm_acc<param::sharp_spherical>(body & particle, std::vector<body *> &nbs) {

    point_t a_apm = 0.0;

    return a_apm;

}
#ifdef apm_type
#define sph_compute_apm_acc compute_apm_acc<param::apm_type>
#else
compute_apm_acc_t sph_compute_apm_acc = nullptr;
#endif
/**
 * @brief      APM selector
 */
void
select() {
  using namespace param;
#ifndef apm_type
  switch(apm_type) {
    case(zero_apm):
      sph_compute_apm_acc = compute_apm_acc<zero_apm>;
      break;
    case(kn_ejecta):
      sph_compute_apm_acc = compute_apm_acc<kn_ejecta>;
      break;
    case(sharp_spherical):
      sph_compute_apm_acc = compute_apm_acc<sharp_spherical>;
      break;
    default:
      log_fatal("Bad APM parameter" << std::endl);
  }
#endif

} // select()

} // namespace apm

