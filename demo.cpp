#include <iostream>

#include "PGA3D.hpp"

int main(int argc, char **argv) {
  static const char *basis[] = {"1",    "e0",   "e1",   "e2",   "e3",  "e01",
                                "e02",  "e03",  "e12",  "e31",  "e23", "e021",
                                "e013", "e032", "e123", "e0123"};

  PGA3D orig = e123;

  // ^ meet
  // & join (first object is the "endpoint")
  // * geometric product

  // Elements of the even subalgebra (scalar + bivector + pss) of unit length
  // are motors PGA3D rot = rotor(PI / 2.0f, e1 * e2);2

  // Elements of the basis of PGA3D
  std::cout << "--- BASIS ELEMENTS ---" << std::endl;
  std::cout << "1:              " << PGA3D() + 1 << std::endl;
  std::cout << "e0:             " << e0 << std::endl;
  std::cout << "x=0 plane = e1: " << e1 << std::endl;
  std::cout << "y=0 plane = e2: " << e2 << std::endl;
  std::cout << "z=0 plane = e3: " << e3 << std::endl;
  std::cout << "e01:            " << e01 << std::endl;
  std::cout << "e02:            " << e02 << std::endl;
  std::cout << "e03:            " << e03 << std::endl;
  std::cout << "z-axis = e12:   " << e12 << std::endl;
  std::cout << "y-axis = e31:   " << e31 << std::endl;
  std::cout << "x-axis = e23:   " << e23 << std::endl;
  std::cout << "e021:           " << e021 << std::endl;
  std::cout << "e013:           " << e013 << std::endl;
  std::cout << "e032:           " << e032 << std::endl;
  std::cout << "origin = e123:  " << e123 << std::endl;
  std::cout << "I = e0123:      " << I << std::endl;
  std::cout << std::endl;

  std::cout << "--- Creating geometric primitives ---" << std::endl;

  // Point from (x,y,z) coordinates
  PGA3D P1 = point(3.0, 4.0, 5.0);
  PGA3D P2 = point(1.0, 2.0, 3.0);
  std::cout << "point P1: " << P1 << std::endl;
  std::cout << "point P2: " << P2 << std::endl;

  // A line through P1 and P2, using the regressive product
  PGA3D l1 = P1 & P2;
  std::cout << "line l: " << l1 << std::endl;

  // The y-direction unit vector from two points
  PGA3D py = point(0.0, 1.0, 0.0);
  PGA3D ax_y_2 = py & orig;
  std::cout << "y axis: " << ax_y_2 << std::endl;

  // We can also easily create points and join them into a line using the
  // regressive (vee, &) product.

  // line in direction
  // PGA3D l1 = point(6.0, 8.0, 14.0) & px;
  // std::cout << "l1: " << l1 << std::endl;

  // plane from normal through point
  PGA3D p3 = point(1, 0, 0);
  PGA3D v1 = p3 & orig;
  PGA3D pn = plane_normal_to_line_through_point(v1, p3);
  std::cout << "pn: " << pn << std::endl;

  // intersect line and plane
  std::cout << "\nLine - plane intersection" << std::endl;
  std::cout << "line (x-axis): " << v1 << std::endl;
  std::cout << "plane (x=0): " << e1 << std::endl;
  PGA3D isec = e1 ^ v1;
  std::cout << "isec (origin): " << isec << std::endl;

  // prepare test for outer product
  // std::cout << "\nPrepare test for outer product" << std::endl;
  // PGA3D mvec1[] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
  // std::cout << mvec1 << std::endl;

  return 0;
}
