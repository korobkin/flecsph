#include "gtest/gtest.h"

#include <cmath>
#include <iostream>
#include <mpi.h>

#include "eforce.h"
#include "tree.h"
#include "params.h"

// default params:
// phi = S (r - R)^n
// S = 1.0E12
// R = 1.0
// n = 5.0
// box = {length=3.0, width=height=1.0}

TEST(eforces, spherical_wall) {

    point_t p{2.0, 0.0, 0.0};
    body b;
    b.set_coordinates(p);

    external_force::force_spherical_wall f;

    double pot = f.potential(b.coordinates());
    point_t acc = f.acceleration(b);

    double exp_pot = 1.0E12;
    double exp_a0 = -5.0E12;
    double exp_a1 = 0.0;
    double exp_a2 = 0.0;

    EXPECT_DOUBLE_EQ(pot, exp_pot);
    EXPECT_DOUBLE_EQ(acc[0], exp_a0);
    EXPECT_DOUBLE_EQ(acc[1], exp_a1);
    EXPECT_DOUBLE_EQ(acc[2], exp_a2);
}

TEST(eforces, square_well)
{
    point_t p{2.0, 0.0, 0.0};

    external_force::force_square_well<0> f;

    double pot = f.potential(b.coordinates());
    point_t acc = f.acceleration(b);

    double exp_pot = 3.125E10;
    double exp_acc = 3.125E11;

    EXPECT_DOUBLE_EQ(pot, exp_pot);
    EXPECT_DOUBLE_EQ(acc, exp_acc)
}

