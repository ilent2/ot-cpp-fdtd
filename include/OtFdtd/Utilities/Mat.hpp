/* Utilities/Mat.hpp - A 3-dimensional matrix type for FDTD simulations.
 *
 * Copyright (C) 2016 Isaac Lenton (aka ilent2)
 */

#ifndef UTILITIES_MAT_HPP
#define UTILITIES_MAT_HPP

#include "Vec.hpp"

namespace Utilities {

template <typename T>
class Mat3t
{
public:

  /** Default constructor: initialize to all zeros.
   *
   * Declaring a separate default constructor rather than default arguments
   * to avoid accidental incomplete construction.
   */
  Mat3t(void) : Mat3t(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0) {}

  /** Construct a matrix from 9 values of type T. */
  Mat3t(T _xx, T _xy, T _xz,
        T _yx, T _yy, T _yz,
        T _zx, T _zy, T _zz)
    : xx(_xx), xy(_xy), xz(_xz),
      yx(_yx), yy(_yy), yz(_yz),
      zx(_zx), zy(_zy), zz(_zz) {}

  T xx, xy, xz, yx, yy, yz, zx, zy, zz;

  /** Return the identity matrix. */
  static Mat3t<T> Eye(void)
  {
    return Mat3t<T>(1.0, 0.0, 0.0,
                    0.0, 1.0, 0.0,
                    0.0, 0.0, 1.0);
  }

  /** Return a matrix of all ones. */
  static Mat3t<T> Ones(void)
  {
    return Mat3t<T>(1.0, 1.0, 1.0,
                    1.0, 1.0, 1.0,
                    1.0, 1.0, 1.0);
  }

  Mat3t<T> inverse(void) const
  {
    double det = xx*(yy*zz - yz*zy) - xy*(yx*zz - yz*zx) + xz*(yx*zy - yy*zx);
    return Mat3t<T>(yy*zz - yz*zy, xz*zy - xy*zz, xy*yz - yy*xz,
                    yz*zx - yx*zz, zz*xx - zx*xz, xz*yx - yz*xx,
                    yx*zy - zx*yy, xy*zx - xx*zy, xx*yy - yx*xy) / det;
  }

  Mat3t<T> t(void) const
  {
    return Mat3t(xx, yx, zx,
                 xy, yy, zy,
                 xz, yz, zz);
  }

  Vec3t<T> row(int r)
  {
    if (r == 2) return Vec3t<T>(xx, xy, xz);
    if (r == 1) return Vec3t<T>(yx, yy, yz);
    if (r == 0) return Vec3t<T>(zx, zy, zz);
    throw std::logic_error("3-D matrix has only 3 dimensions");
  }

  Vec3t<T> col(int r)
  {
    if (r == 2) return Vec3t<T>(xx, yx, zx);
    if (r == 1) return Vec3t<T>(xy, yy, zy);
    if (r == 0) return Vec3t<T>(xz, yz, zz);
    throw std::logic_error("3-D matrix has only 3 dimensions");
  }

  Mat3t<T> cross(const Utilities::Vec3t<T>& vec)
  {
    return (*this) * Mat3t<T>(0.0, vec.z, -vec.y,
                              -vec.z, 0.0, vec.x,
                              vec.y, -vec.x, 0.0);
  }

  Mat3t<T> operator+ (const Mat3t<T>& o) const
  {
    return Mat3t<T>(xx + o.xx, xy + o.xy, xz + o.xz,
                    yx + o.yx, yy + o.yy, yz + o.yz,
                    zx + o.zx, zy + o.zy, zz + o.zz);
  }

  Mat3t<T> operator- (const Mat3t<T>& o) const
  {
    return Mat3t<T>(xx - o.xx, xy - o.xy, xz - o.xz,
                    yx - o.yx, yy - o.yy, yz - o.yz,
                    zx - o.zx, zy - o.zy, zz - o.zz);
  }
};

template <typename T>
Vec3t<T> operator* (Vec3t<T> a, Mat3t<T> b)
{
  return Vec3t<T>(a.x*b.xx + a.y*b.yx + a.z*b.zx,
                  a.x*b.xy + a.y*b.yy + a.z*b.zy,
                  a.x*b.xz + a.y*b.yz + a.z*b.zz);
}

template <typename T>
Vec3t<T> operator* (Mat3t<T> a, Vec3t<T> b)
{
  return Vec3t<T>(a.xx*b.x + a.xy*b.y + a.xz*b.z,
                  a.yx*b.x + a.yy*b.y + a.yz*b.z,
                  a.zx*b.x + a.zy*b.y + a.zz*b.z);
}

template <typename T>
Mat3t<T> operator* (Mat3t<T> a, Mat3t<T> b)
{
  return Mat3t<T>(a.xx*b.xx + a.xy*b.yx + a.xz*b.zx,
                    a.xx*b.xy + a.xy*b.yy + a.xz*b.zy,
                    a.xx*b.xz + a.xy*b.yz + a.xz*b.zz,
                  a.yx*b.xx + a.yy*b.yx + a.yz*b.zx,
                    a.yx*b.xy + a.yy*b.yy + a.yz*b.zy,
                    a.yx*b.xz + a.yy*b.yz + a.yz*b.zz,
                  a.zx*b.xx + a.zy*b.yx + a.zz*b.zx,
                    a.zx*b.xy + a.zy*b.yy + a.zz*b.zy,
                    a.zx*b.xz + a.zy*b.yz + a.zz*b.zz);
}

template <typename T>
Mat3t<T> operator/ (Mat3t<T> a, T b)
{
  return Mat3t<T>(a.xx/b, a.xy/b, a.xz/b,
                  a.yx/b, a.yy/b, a.yz/b,
                  a.zx/b, a.zy/b, a.zz/b);
}

template <class O, class T>
O& operator<<(O& ostream, const Mat3t<T>& a)
{
  ostream << a.xx << '\t' << a.xy << '\t' << a.xz << '\n';
  ostream << a.yx << '\t' << a.yy << '\t' << a.yz << '\n';
  ostream << a.zx << '\t' << a.zy << '\t' << a.zz << '\n';
  return ostream;
}


typedef Mat3t<double> Mat3d;

} // namespace Utilities

#endif // #ifndef UTILITIES_MAT_HPP

