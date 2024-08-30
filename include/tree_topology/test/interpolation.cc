#include "gtest/gtest.h"

#include <iomanip>
#include <iostream>

#include "interp.h"

namespace flecsi {
namespace execution {
void
driver(int, char **) {}
} // namespace execution
} // namespace flecsi

TEST(interpolation, 1d_linear) {
  double ys[] = {0.0, 0.1, 0.3, 0.6, 1.0};
  interp::linear_interpolator_1d f{0, 1, ys, sizeof(ys)/sizeof(double)};
  //std::cout << "f(0.215) = " << f(0.215) << std::endl;
  EXPECT_DOUBLE_EQ(f(0.215), 0.086);
  EXPECT_DOUBLE_EQ(f(0.), 0.);
  EXPECT_DOUBLE_EQ(f(1.), 1.);
  EXPECT_DOUBLE_EQ(f(-1.),0.);
  EXPECT_DOUBLE_EQ(f(10.),1.);
}

TEST(interpolation, 1d_cubic) {
  double ys[] = {0.0, 0.008, 0.064, 0.216, 0.512, 1.0};
  interp::cubic_interpolator_1d f{0, 1, ys, sizeof(ys)/sizeof(double)};
  //std::cout << "f(0.215) = " << f(0.215) << std::endl;
  EXPECT_NEAR(f(0.215), 0.215*0.215*0.215, 1e-12);
  EXPECT_DOUBLE_EQ(f(2.0), 1.0);
  EXPECT_DOUBLE_EQ(f(-2.0), 0.0);
  EXPECT_DOUBLE_EQ(f(0.), 0.);
  EXPECT_DOUBLE_EQ(f(1.), 1.);
}

TEST(interpolation, 1d_linear_nonuniform_grid) {
  double xs[] = {0.0, 0.5, 0.8, 0.9, 1.0};
  double ys[] = {0.0, 0.1, 0.3, 0.6, 1.0};
  interp::linear_interpolator_1d_nug f{xs, ys, sizeof(ys)/sizeof(double)};
  //std::cout << "f(0.215) = " << f(0.215) << std::endl;
  EXPECT_DOUBLE_EQ(f(0.215), 0.043);
  EXPECT_DOUBLE_EQ(f(2.0), 1.0);
  EXPECT_DOUBLE_EQ(f(-2.0), 0.0);
  EXPECT_DOUBLE_EQ(f(0.), 0.);
  EXPECT_DOUBLE_EQ(f(1.), 1.);
}

TEST(interpolation, 1d_cubic_nonuniform_grid) {
  double xs[] = {0.0, 0.5, 0.7, 0.9, 1.0};
  double ys[] = {0.0, 0.125, 0.343, 0.729, 1.0};
  interp::cubic_interpolator_1d_nug f{xs, ys, sizeof(ys)/sizeof(double)};
  //std::cout << "f(0.215) = " << f(0.215) << std::endl;
  //for (double x = 0; x < 1; x += 0.01) printf("%5.2f %10.7f\n", x, f(x));
  EXPECT_NEAR(f(0.215), 0.215*0.215*0.215, 1e-12);
  EXPECT_DOUBLE_EQ(f(2.0), 1.0);
  EXPECT_DOUBLE_EQ(f(-2.0), 0.0);
  EXPECT_DOUBLE_EQ(f(0.), 0.);
  EXPECT_DOUBLE_EQ(f(1.), 1.);
}

TEST(interpolation, 1d_cubic_nonuniform_grid_float) {
  using namespace interp;
  double xs[] = {0.0, 0.5, 0.7, 0.9, 1.0};
  float ys[] = {0.0, 0.125, 0.343, 0.729, 1.0};
  interpolating_function_u<float, 1, nonuniform_grid, 3>
      f{xs, ys, sizeof(ys)/sizeof(float)};
  //std::cout << "f(0.215) = " << f(0.215) << std::endl;
  //for (double x = 0; x < 1; x += 0.01) printf("%5.2f %10.7f\n", x, f(x));
  EXPECT_NEAR(f(0.215), 0.215*0.215*0.215, 1e-8);
  EXPECT_DOUBLE_EQ(f(2.0), 1.0);
  EXPECT_DOUBLE_EQ(f(-2.0), 0.0);
  EXPECT_DOUBLE_EQ(f(0.), 0.);
  EXPECT_DOUBLE_EQ(f(1.), 1.);
}

TEST(interpolation, 2d_linear) {
  double ys[] = {0.1, -0.1, -0.3,
                 0.6,  1.0,  2.4};
  interp::linear_interpolator_2d f{-1., 1., 0., 1., ys, 3, 2};
  //for (double x = -1; x <= 1. + .001; x += 0.01) {
  //  for (double y = 0; y <= 1. + .001; y += 0.01)
  //    printf("%5.2f %5.2f %12.9f\n", x, y, f(x, y));
  //  printf("\n");
  //}
  //std::cout << "f(0.215) = " << f(0.215) << std::endl;
  EXPECT_DOUBLE_EQ(f(-0.4, 0.215), 0.1649);
  EXPECT_DOUBLE_EQ(f(-1., 0.), 0.1);
  EXPECT_DOUBLE_EQ(f(-1., 1.), 0.6);
  EXPECT_DOUBLE_EQ(f( 1., 0.),-0.3);
  EXPECT_DOUBLE_EQ(f( 1., 1.), 2.4);
}

TEST(interpolation, 2d_linear_nonuniform) {
  double xs[] = {-1.0, 0.5, 1.0}; // x-grid
  double ys[] = { 0.0, 0.1, 0.9, 1.0}; // y-grid
  double fs[] = {0.0, 0.1, 0.3, // column-major order
                 0.1,-1.0, 2.1,
                 0.3,-1.0, 2.0,
                 0.6, 0.3, 2.0};
  interp::linear_interpolator_2d_nug f{xs, ys, fs, 3, 4};
  //std::cout << "f(-0.1, 0.215) = " << f(-0.1, 0.215) << std::endl;
  for (double x = -1; x <= 1. + .001; x += 0.01) {
    for (double y = 0; y <= 1. + .001; y += 0.01)
      printf("%5.2f %5.2f %12.9f\n", x, y, f(x, y));
    printf("\n");
  }
  //std::cout << "f(0.215) = " << f(0.215) << std::endl;
}

