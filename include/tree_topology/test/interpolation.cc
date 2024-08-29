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

