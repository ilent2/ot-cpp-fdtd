/* VecRtp.hpp - Vector types for spherical coordinates.
 *
 * Copyright (C) 2016 Isaac Lenton (aka ilent2)
 */

#ifndef UTILITIES_VEC_RTP_HPP
#define UTILITIES_VEC_RTP_HPP

#include <complex>
#include "Vec.hpp"        // For type conversion, is there a better way?

namespace Utilities {

/** A spherical coordinate type. */
template <typename T>
class VecRtpt
{
public:

  /** Create a new vector. */
  VecRtpt(T _r=0.0, T _t=0.0, T _p=0.0) : r(_r), t(_t), p(_p) {}

  /** Create a new vector from a Vec3t<T>. */
  VecRtpt(const Vec3t<T>& other)
  {
    double rxy = pow(pow(other.x, 2.0) + pow(other.y, 2.0), 0.5);

    r = other.length();
    t = atan2(rxy, other.z);                              // already positive
    p = fmod(atan2(other.y, other.x)+2.0*M_PI, 2.0*M_PI); // make positive
  }

  /** Convert to a Vec3t<T>. */
  operator Vec3t<T>() const
  {
    return Vec3t<T>(r*cos(p)*sin(t), r*sin(p)*sin(t), r*cos(t));
  }

  /** Return the length of the vector. */
  T length(void) const { return r; }

  /** Compute the length^2 of the vector. */
  T length2(void) const { return r*r; }

  /** Return a normalized copy of ourselves. */
  VecRtpt<T> normalized(void) const
  {
    return VecRtpt<T>(1.0, t, p);
  }

  /** Add things component wise, useful for adding derivative values. */
  VecRtpt<T> add_components(const VecRtpt<T>& o) const
  {
    return VecRtpt<T>(r + o.r, t + o.t, p + o.p);
  }

  /** Multiply things component wise, useful for derivative values. */
  template <typename O>
  VecRtpt<T> mul_components(const O& o) const
  {
    return VecRtpt<T>(r * o, t * o, p * o);
  }

  /** Multiply things component wise, useful for derivative values. */
  VecRtpt<T> mul_components(const VecRtpt<T>& o) const
  {
    return VecRtpt<T>(r * o.r, t * o.t, p * o.p);
  }

  T r;        ///< The length of the vector (radial)
  T t;        ///< The angle to the z axis (theta)
  T p;        ///< The angle around the z-axis from the x-axis (phi)
};

typedef VecRtpt<double> VecRtp;
typedef VecRtpt<std::complex<double>> VecRtpc;

} // namespace Utilities

#endif // #ifndef UTILITIES_VEC_RTP_HPP

