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

TEST(interpolation, sanity) {

  using namespace std;
  using namespace interp;

  cout << "--- Testing 1D interpolation: ---" << endl;
  double ys[5] = {0.0, 0.1, 0.3, 0.6, 1.0};
  interp::interpolating_function_1d f{0, 1, ys, sizeof(ys)/sizeof(double)};
  std::cout << "f(0.215) = " << f(0.215) << std::endl; 
  std::cout << "f(0.0) = " << f(0.0) << std::endl; 
  std::cout << "f(1.0) = " << f(1.0) << std::endl; 
  std::cout << "f(-1.0) = " << f(-1.0) << std::endl; 
  std::cout << "f(10.0) = " << f(10.0) << std::endl; 
  EXPECT_DOUBLE_EQ(f(0.215), 0.086);
  EXPECT_DOUBLE_EQ(f(0.), 0.);
  EXPECT_DOUBLE_EQ(f(1.), 1.);
  EXPECT_DOUBLE_EQ(f(-1.),0.);
  EXPECT_DOUBLE_EQ(f(10.),1.);
}

