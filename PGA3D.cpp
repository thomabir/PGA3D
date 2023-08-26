#include "PGA3D.hpp"

#include <array>
#include <cmath>
#include <iostream>

// Default constructor: zero vector
PGA3D::PGA3D() { mvec.fill(0.0); }

// Constructor: Exactly one index non-zero
PGA3D::PGA3D(double f, int idx = 0) {
  mvec.fill(0.0);
  mvec[idx] = f;
}

std::array<double, 16> PGA3D::get_mvec() { return mvec; }

// Return a single element of vector
double &PGA3D::operator[](size_t idx) { return mvec[idx]; }
const double &PGA3D::operator[](size_t idx) const { return mvec[idx]; }

// Print mvec
std::ostream &operator<<(std::ostream &out, const PGA3D &v) {
  out << '(';
  for (int i = 0; i < 15; i++) {
    out << v.mvec[i] << ", ";
  }

  out << v.mvec[15] << ')';
  return out;
}

// Unary minus
PGA3D operator-(const PGA3D &a) {
  PGA3D res;
  for (int i = 0; i < 16; i++) {
    res[i] = -a[i];
  }
  return res;
}

// Basic linear algebra operations: addition, subtraction, scalar operations

// Vector addition
PGA3D operator+(const PGA3D &a, const PGA3D &b) {
  PGA3D res;
  for (int i = 0; i < 16; i++) {
    res[i] = a[i] + b[i];
  }
  return res;
}

// Vector subtraction
PGA3D operator-(const PGA3D &a, const PGA3D &b) { return a + (-b); }

// Scalar multiplication (from left and from right)
PGA3D operator*(const double &a, const PGA3D &b) {
  PGA3D res;
  for (int i = 0; i < 16; i++) {
    res[i] = a * b[i];
  }
  return res;
}

PGA3D operator*(const PGA3D &a, const double &b) { return b * a; }

// Scalar addition (from left and from right)
PGA3D operator+(const double &a, const PGA3D &b) {
  PGA3D res = b;
  res[0] += a;
  return res;
}

PGA3D operator+(const PGA3D &a, const double &b) { return b + a; }

// Scalar subtraction (from left and from right)
PGA3D operator-(const PGA3D &a, const double &b) { return a + (-b); }

PGA3D operator-(const double &a, const PGA3D &b) { return a + (-b); }

// Projective geometric algebra operations: geometric product, inner product ...

// Reverse the order of the basis blades
PGA3D operator~(const PGA3D &a) {
  PGA3D res;
  res[0] = a[0];
  res[1] = a[1];
  res[2] = a[2];
  res[3] = a[3];
  res[4] = a[4];
  res[5] = -a[5];
  res[6] = -a[6];
  res[7] = -a[7];
  res[8] = -a[8];
  res[9] = -a[9];
  res[10] = -a[10];
  res[11] = -a[11];
  res[12] = -a[12];
  res[13] = -a[13];
  res[14] = -a[14];
  res[15] = a[15];
  return res;
}

// Poincare duality operator
PGA3D operator!(const PGA3D &a) {
  PGA3D res;
  res[0] = a[15];
  res[1] = a[14];
  res[2] = a[13];
  res[3] = a[12];
  res[4] = a[11];
  res[5] = a[10];
  res[6] = a[9];
  res[7] = a[8];
  res[8] = a[7];
  res[9] = a[6];
  res[10] = a[5];
  res[11] = a[4];
  res[12] = a[3];
  res[13] = a[2];
  res[14] = a[1];
  res[15] = a[0];
  return res;
}

// Clifford conjugation
PGA3D PGA3D::Conjugate() {
  PGA3D res;
  res[0] = this->mvec[0];
  res[1] = -this->mvec[1];
  res[2] = -this->mvec[2];
  res[3] = -this->mvec[3];
  res[4] = -this->mvec[4];
  res[5] = -this->mvec[5];
  res[6] = -this->mvec[6];
  res[7] = -this->mvec[7];
  res[8] = -this->mvec[8];
  res[9] = -this->mvec[9];
  res[10] = -this->mvec[10];
  res[11] = this->mvec[11];
  res[12] = this->mvec[12];
  res[13] = this->mvec[13];
  res[14] = this->mvec[14];
  res[15] = this->mvec[15];
  return res;
}

// Main involution
PGA3D PGA3D::Involute() {
  PGA3D res;
  res[0] = this->mvec[0];
  res[1] = -this->mvec[1];
  res[2] = -this->mvec[2];
  res[3] = -this->mvec[3];
  res[4] = -this->mvec[4];
  res[5] = this->mvec[5];
  res[6] = this->mvec[6];
  res[7] = this->mvec[7];
  res[8] = this->mvec[8];
  res[9] = this->mvec[9];
  res[10] = this->mvec[10];
  res[11] = -this->mvec[11];
  res[12] = -this->mvec[12];
  res[13] = -this->mvec[13];
  res[14] = -this->mvec[14];
  res[15] = this->mvec[15];
  return res;
}

// Geometric product
PGA3D operator*(const PGA3D &a, const PGA3D &b) {
  PGA3D res;
  res[0] = b[0] * a[0] + b[2] * a[2] + b[3] * a[3] + b[4] * a[4] - b[8] * a[8] -
           b[9] * a[9] - b[10] * a[10] - b[14] * a[14];
  res[1] = b[1] * a[0] + b[0] * a[1] - b[5] * a[2] - b[6] * a[3] - b[7] * a[4] +
           b[2] * a[5] + b[3] * a[6] + b[4] * a[7] + b[11] * a[8] +
           b[12] * a[9] + b[13] * a[10] + b[8] * a[11] + b[9] * a[12] +
           b[10] * a[13] + b[15] * a[14] - b[14] * a[15];
  res[2] = b[2] * a[0] + b[0] * a[2] - b[8] * a[3] + b[9] * a[4] + b[3] * a[8] -
           b[4] * a[9] - b[14] * a[10] - b[10] * a[14];
  res[3] = b[3] * a[0] + b[8] * a[2] + b[0] * a[3] - b[10] * a[4] -
           b[2] * a[8] - b[14] * a[9] + b[4] * a[10] - b[9] * a[14];
  res[4] = b[4] * a[0] - b[9] * a[2] + b[10] * a[3] + b[0] * a[4] -
           b[14] * a[8] + b[2] * a[9] - b[3] * a[10] - b[8] * a[14];
  res[5] = b[5] * a[0] + b[2] * a[1] - b[1] * a[2] - b[11] * a[3] +
           b[12] * a[4] + b[0] * a[5] - b[8] * a[6] + b[9] * a[7] +
           b[6] * a[8] - b[7] * a[9] - b[15] * a[10] - b[3] * a[11] +
           b[4] * a[12] + b[14] * a[13] - b[13] * a[14] - b[10] * a[15];
  res[6] = b[6] * a[0] + b[3] * a[1] + b[11] * a[2] - b[1] * a[3] -
           b[13] * a[4] + b[8] * a[5] + b[0] * a[6] - b[10] * a[7] -
           b[5] * a[8] - b[15] * a[9] + b[7] * a[10] + b[2] * a[11] +
           b[14] * a[12] - b[4] * a[13] - b[12] * a[14] - b[9] * a[15];
  res[7] = b[7] * a[0] + b[4] * a[1] - b[12] * a[2] + b[13] * a[3] -
           b[1] * a[4] - b[9] * a[5] + b[10] * a[6] + b[0] * a[7] -
           b[15] * a[8] + b[5] * a[9] - b[6] * a[10] + b[14] * a[11] -
           b[2] * a[12] + b[3] * a[13] - b[11] * a[14] - b[8] * a[15];
  res[8] = b[8] * a[0] + b[3] * a[2] - b[2] * a[3] + b[14] * a[4] +
           b[0] * a[8] + b[10] * a[9] - b[9] * a[10] + b[4] * a[14];
  res[9] = b[9] * a[0] - b[4] * a[2] + b[14] * a[3] + b[2] * a[4] -
           b[10] * a[8] + b[0] * a[9] + b[8] * a[10] + b[3] * a[14];
  res[10] = b[10] * a[0] + b[14] * a[2] + b[4] * a[3] - b[3] * a[4] +
            b[9] * a[8] - b[8] * a[9] + b[0] * a[10] + b[2] * a[14];
  res[11] = b[11] * a[0] - b[8] * a[1] + b[6] * a[2] - b[5] * a[3] +
            b[15] * a[4] - b[3] * a[5] + b[2] * a[6] - b[14] * a[7] -
            b[1] * a[8] + b[13] * a[9] - b[12] * a[10] + b[0] * a[11] +
            b[10] * a[12] - b[9] * a[13] + b[7] * a[14] - b[4] * a[15];
  res[12] = b[12] * a[0] - b[9] * a[1] - b[7] * a[2] + b[15] * a[3] +
            b[5] * a[4] + b[4] * a[5] - b[14] * a[6] - b[2] * a[7] -
            b[13] * a[8] - b[1] * a[9] + b[11] * a[10] - b[10] * a[11] +
            b[0] * a[12] + b[8] * a[13] + b[6] * a[14] - b[3] * a[15];
  res[13] = b[13] * a[0] - b[10] * a[1] + b[15] * a[2] + b[7] * a[3] -
            b[6] * a[4] - b[14] * a[5] - b[4] * a[6] + b[3] * a[7] +
            b[12] * a[8] - b[11] * a[9] - b[1] * a[10] + b[9] * a[11] -
            b[8] * a[12] + b[0] * a[13] + b[5] * a[14] - b[2] * a[15];
  res[14] = b[14] * a[0] + b[10] * a[2] + b[9] * a[3] + b[8] * a[4] +
            b[4] * a[8] + b[3] * a[9] + b[2] * a[10] + b[0] * a[14];
  res[15] = b[15] * a[0] + b[14] * a[1] + b[13] * a[2] + b[12] * a[3] +
            b[11] * a[4] + b[10] * a[5] + b[9] * a[6] + b[8] * a[7] +
            b[7] * a[8] + b[6] * a[9] + b[5] * a[10] - b[4] * a[11] -
            b[3] * a[12] - b[2] * a[13] - b[1] * a[14] + b[0] * a[15];
  return res;
}

// Outer product, meet, wedge, exterior product, intersection
PGA3D operator^(const PGA3D &a, const PGA3D &b) {
  PGA3D res;
  res[0] = b[0] * a[0];
  res[1] = b[1] * a[0] + b[0] * a[1];
  res[2] = b[2] * a[0] + b[0] * a[2];
  res[3] = b[3] * a[0] + b[0] * a[3];
  res[4] = b[4] * a[0] + b[0] * a[4];
  res[5] = b[5] * a[0] + b[2] * a[1] - b[1] * a[2] + b[0] * a[5];
  res[6] = b[6] * a[0] + b[3] * a[1] - b[1] * a[3] + b[0] * a[6];
  res[7] = b[7] * a[0] + b[4] * a[1] - b[1] * a[4] + b[0] * a[7];
  res[8] = b[8] * a[0] + b[3] * a[2] - b[2] * a[3] + b[0] * a[8];
  res[9] = b[9] * a[0] - b[4] * a[2] + b[2] * a[4] + b[0] * a[9];
  res[10] = b[10] * a[0] + b[4] * a[3] - b[3] * a[4] + b[0] * a[10];
  res[11] = b[11] * a[0] - b[8] * a[1] + b[6] * a[2] - b[5] * a[3] -
            b[3] * a[5] + b[2] * a[6] - b[1] * a[8] + b[0] * a[11];
  res[12] = b[12] * a[0] - b[9] * a[1] - b[7] * a[2] + b[5] * a[4] +
            b[4] * a[5] - b[2] * a[7] - b[1] * a[9] + b[0] * a[12];
  res[13] = b[13] * a[0] - b[10] * a[1] + b[7] * a[3] - b[6] * a[4] -
            b[4] * a[6] + b[3] * a[7] - b[1] * a[10] + b[0] * a[13];
  res[14] = b[14] * a[0] + b[10] * a[2] + b[9] * a[3] + b[8] * a[4] +
            b[4] * a[8] + b[3] * a[9] + b[2] * a[10] + b[0] * a[14];
  res[15] = b[15] * a[0] + b[14] * a[1] + b[13] * a[2] + b[12] * a[3] +
            b[11] * a[4] + b[10] * a[5] + b[9] * a[6] + b[8] * a[7] +
            b[7] * a[8] + b[6] * a[9] + b[5] * a[10] - b[4] * a[11] -
            b[3] * a[12] - b[2] * a[13] - b[1] * a[14] + b[0] * a[15];
  return res;
}

// Regressive product, join, vee
PGA3D operator&(const PGA3D &a, const PGA3D &b) {
  PGA3D res;
  res[15] = 1 * (a[15] * b[15]);
  res[14] = -1 * (a[14] * -1 * b[15] + a[15] * b[14] * -1);
  res[13] = -1 * (a[13] * -1 * b[15] + a[15] * b[13] * -1);
  res[12] = -1 * (a[12] * -1 * b[15] + a[15] * b[12] * -1);
  res[11] = -1 * (a[11] * -1 * b[15] + a[15] * b[11] * -1);
  res[10] = 1 * (a[10] * b[15] + a[13] * -1 * b[14] * -1 -
                 a[14] * -1 * b[13] * -1 + a[15] * b[10]);
  res[9] = 1 * (a[9] * b[15] + a[12] * -1 * b[14] * -1 -
                a[14] * -1 * b[12] * -1 + a[15] * b[9]);
  res[8] = 1 * (a[8] * b[15] + a[11] * -1 * b[14] * -1 -
                a[14] * -1 * b[11] * -1 + a[15] * b[8]);
  res[7] = 1 * (a[7] * b[15] + a[12] * -1 * b[13] * -1 -
                a[13] * -1 * b[12] * -1 + a[15] * b[7]);
  res[6] = 1 * (a[6] * b[15] - a[11] * -1 * b[13] * -1 +
                a[13] * -1 * b[11] * -1 + a[15] * b[6]);
  res[5] = 1 * (a[5] * b[15] + a[11] * -1 * b[12] * -1 -
                a[12] * -1 * b[11] * -1 + a[15] * b[5]);
  res[4] = 1 * (a[4] * b[15] - a[7] * b[14] * -1 + a[9] * b[13] * -1 -
                a[10] * b[12] * -1 - a[12] * -1 * b[10] + a[13] * -1 * b[9] -
                a[14] * -1 * b[7] + a[15] * b[4]);
  res[3] = 1 * (a[3] * b[15] - a[6] * b[14] * -1 - a[8] * b[13] * -1 +
                a[10] * b[11] * -1 + a[11] * -1 * b[10] - a[13] * -1 * b[8] -
                a[14] * -1 * b[6] + a[15] * b[3]);
  res[2] = 1 * (a[2] * b[15] - a[5] * b[14] * -1 + a[8] * b[12] * -1 -
                a[9] * b[11] * -1 - a[11] * -1 * b[9] + a[12] * -1 * b[8] -
                a[14] * -1 * b[5] + a[15] * b[2]);
  res[1] = 1 * (a[1] * b[15] + a[5] * b[13] * -1 + a[6] * b[12] * -1 +
                a[7] * b[11] * -1 + a[11] * -1 * b[7] + a[12] * -1 * b[6] +
                a[13] * -1 * b[5] + a[15] * b[1]);
  res[0] = 1 * (a[0] * b[15] + a[1] * b[14] * -1 + a[2] * b[13] * -1 +
                a[3] * b[12] * -1 + a[4] * b[11] * -1 + a[5] * b[10] +
                a[6] * b[9] + a[7] * b[8] + a[8] * b[7] + a[9] * b[6] +
                a[10] * b[5] - a[11] * -1 * b[4] - a[12] * -1 * b[3] -
                a[13] * -1 * b[2] - a[14] * -1 * b[1] + a[15] * b[0]);
  return res;
}

// Inner product, dot product
PGA3D operator|(const PGA3D &a, const PGA3D &b) {
  PGA3D res;
  res[0] = b[0] * a[0] + b[2] * a[2] + b[3] * a[3] + b[4] * a[4] - b[8] * a[8] -
           b[9] * a[9] - b[10] * a[10] - b[14] * a[14];
  res[1] = b[1] * a[0] + b[0] * a[1] - b[5] * a[2] - b[6] * a[3] - b[7] * a[4] +
           b[2] * a[5] + b[3] * a[6] + b[4] * a[7] + b[11] * a[8] +
           b[12] * a[9] + b[13] * a[10] + b[8] * a[11] + b[9] * a[12] +
           b[10] * a[13] + b[15] * a[14] - b[14] * a[15];
  res[2] = b[2] * a[0] + b[0] * a[2] - b[8] * a[3] + b[9] * a[4] + b[3] * a[8] -
           b[4] * a[9] - b[14] * a[10] - b[10] * a[14];
  res[3] = b[3] * a[0] + b[8] * a[2] + b[0] * a[3] - b[10] * a[4] -
           b[2] * a[8] - b[14] * a[9] + b[4] * a[10] - b[9] * a[14];
  res[4] = b[4] * a[0] - b[9] * a[2] + b[10] * a[3] + b[0] * a[4] -
           b[14] * a[8] + b[2] * a[9] - b[3] * a[10] - b[8] * a[14];
  res[5] = b[5] * a[0] - b[11] * a[3] + b[12] * a[4] + b[0] * a[5] -
           b[15] * a[10] - b[3] * a[11] + b[4] * a[12] - b[10] * a[15];
  res[6] = b[6] * a[0] + b[11] * a[2] - b[13] * a[4] + b[0] * a[6] -
           b[15] * a[9] + b[2] * a[11] - b[4] * a[13] - b[9] * a[15];
  res[7] = b[7] * a[0] - b[12] * a[2] + b[13] * a[3] + b[0] * a[7] -
           b[15] * a[8] - b[2] * a[12] + b[3] * a[13] - b[8] * a[15];
  res[8] = b[8] * a[0] + b[14] * a[4] + b[0] * a[8] + b[4] * a[14];
  res[9] = b[9] * a[0] + b[14] * a[3] + b[0] * a[9] + b[3] * a[14];
  res[10] = b[10] * a[0] + b[14] * a[2] + b[0] * a[10] + b[2] * a[14];
  res[11] = b[11] * a[0] + b[15] * a[4] + b[0] * a[11] - b[4] * a[15];
  res[12] = b[12] * a[0] + b[15] * a[3] + b[0] * a[12] - b[3] * a[15];
  res[13] = b[13] * a[0] + b[15] * a[2] + b[0] * a[13] - b[2] * a[15];
  res[14] = b[14] * a[0] + b[0] * a[14];
  res[15] = b[15] * a[0] + b[0] * a[15];
  return res;
}

double PGA3D::norm() { return sqrt(std::abs(((*this) * Conjugate()).mvec[0])); }
double PGA3D::inorm() { return (!(*this)).norm(); }
PGA3D PGA3D::normalized() { return (*this) * (1 / norm()); }

// basis vectors

// vectors (= "plane ingredients")
const PGA3D e0(1.0f, 1);  // ehh? See plane_from_homogeneous for enlightenment
const PGA3D e1(1.0f, 2);  // x=0 plane
const PGA3D e2(1.0f, 3);  // y=0 plane
const PGA3D e3(1.0f, 4);  // z=0 plane

// trivectors (= "point ingredients")
const PGA3D e123 =
    e1 ^ e2 ^ e3;  // point at origin = intersection of coordinate planes
const PGA3D e032 = e0 ^ e3 ^ e2;  // unit direction vector in the x direction
const PGA3D e013 = e0 ^ e1 ^ e3;  //                              y
const PGA3D e021 = e0 ^ e2 ^ e1;  //                              z

// bivectors (= lines), which are the duals of vectors (?)
const PGA3D e12 = e1 ^ e2;  // z-axis
const PGA3D e31 = e3 ^ e1;  // y-axis
const PGA3D e23 = e2 ^ e3;  // x-axis
const PGA3D e01 = e0 ^ e1;  // ?
const PGA3D e02 = e0 ^ e2;  // ??
const PGA3D e03 = e0 ^ e3;  // ???

// pseudo scalar
const PGA3D I = e0 ^ e1 ^ e2 ^ e3;  // like sqrt(-1) in complex numbers

// motors

// rotor (Euclidean line)
PGA3D rotor(double angle, PGA3D &line) {
  return cos(angle / 2.0) + sin(angle / 2.0) * line.normalized();
}

// translator (Ideal line)
PGA3D translator(double dist, const PGA3D &line) {
  return 1.0 + dist / 2.0 * line;
}

// ways to make a point

// from x,y,z coordinates
PGA3D point(double x, double y, double z) {
  return e123 + x * e032 + y * e013 + z * e021;
}

// ways to make a line
PGA3D lineFromPoints(const PGA3D &pt1, const PGA3D &pt2) { return pt1 & pt2; }
PGA3D lineFromPlanes(const PGA3D &pn1, const PGA3D &pn2) { return pn1 ^ pn2; }
PGA3D lineThroughPointPerpendicularToPlane(const PGA3D &pt, const PGA3D &pn) {
  return pt | pn;
}
PGA3D lineThroughPointParallelToLine(const PGA3D &pt, const PGA3D &ln) {
  return pt | ln * pt;
}
PGA3D lineThroughPointPerpendicularToLine(const PGA3D &pt, const PGA3D &ln) {
  return (pt | ln * pt) & pt;
}

// ways to make a plane

// plane from homogenous equation ax + by + cz + d = 0
PGA3D plane_from_homogeneous(double a, double b, double c, double d) {
  return a * e1 + b * e2 + c * e3 + d * e0;
}

PGA3D plane_normal_to_line_through_point(const PGA3D &l, const PGA3D &P) {
  return l | P;
}

// distance

// angle

PGA3D direction(double x, double y, double z) {
  return x * e032 + y * e013 + z * e021;
}

PGA3D rotate(const PGA3D &object, const PGA3D &rotor) {
  return rotor * object * ~rotor;
}

// generating other geometric shapes
PGA3D circle(double t, double radius, PGA3D line) {
  return rotor(t * 2.0f * M_PI, line) * translator(radius, e1 * e0);
}

PGA3D torus(double s, double t, double r1, PGA3D l1, double r2, PGA3D l2) {
  return circle(s, r2, l2) * circle(t, r1, l1);
}

// PGA3D point_on_torus(double s, double t)
// {
//   PGA3D to = torus(s, t, 0.25f, e1 * e2, 0.6f, e1 * e3);
//   return to * e123 * ~to;
// }
