#pragma once

#include <stddef.h>  // for size_t

#include <array>
#include <iosfwd>  // for ostream

class PGA3D {
 private:
  std::array<double, 16> mvec;

 public:
  PGA3D();
  PGA3D(double f, int idx);

  double &operator[](size_t idx);
  const double &operator[](size_t idx) const;

  std::array<double, 16> get_mvec();
  std::array<double, 4> get_4point();

  PGA3D Conjugate();
  PGA3D Involute();

  double norm();
  double inorm();
  PGA3D normalized();

  friend std::ostream &operator<<(std::ostream &out, const PGA3D &v);
  friend bool operator==(const PGA3D &a, const PGA3D &b);
};

// algebra
PGA3D operator~(const PGA3D &a);
PGA3D operator!(const PGA3D &a);

PGA3D operator*(const PGA3D &a, const PGA3D &b);
PGA3D operator^(const PGA3D &a, const PGA3D &b);
PGA3D operator&(const PGA3D &a, const PGA3D &b);
PGA3D operator|(const PGA3D &a, const PGA3D &b);
PGA3D operator+(const PGA3D &a, const PGA3D &b);
PGA3D operator-(const PGA3D &a, const PGA3D &b);

PGA3D operator*(const double &a, const PGA3D &b);
PGA3D operator*(const PGA3D &a, const double &b);
PGA3D operator+(const double &a, const PGA3D &b);
PGA3D operator+(const PGA3D &a, const double &b);
PGA3D operator-(const PGA3D &a, const double &b);

bool operator==(const PGA3D &a, const PGA3D &b);
bool operator!=(const PGA3D &a, const PGA3D &b);

// vectors (planes)
extern const PGA3D e0, e1, e2, e3;

// trivectors (= points)
extern const PGA3D e123, e032, e013, e021;

// bivectors (= lines)
extern const PGA3D e12, e31, e23;
extern const PGA3D e01, e02, e03;

// pseudo scalar
extern const PGA3D I;

// helpers
PGA3D rotor(double angle, PGA3D &line);
PGA3D translator(double dist, PGA3D &line);
PGA3D plane_from_homogeneous(double a, double b, double c, double d);
PGA3D point(double x, double y, double z);
PGA3D lineFromPoints(const PGA3D &p1, const PGA3D &p2);
PGA3D direction(double x, double y, double z);
PGA3D plane_normal_to_line_through_point(const PGA3D &l, const PGA3D &P);
PGA3D rotate(const PGA3D &object, const PGA3D &rotor);
