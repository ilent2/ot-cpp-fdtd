/** EmField.hpp - Common 1, 2 and 3 dimensional EM field elements.
 *
 * Copyright (C) 2016 Isaac Lenton (aka ilent2).
 */

#ifndef UTILITIES_VEC_HPP
#define UTILITIES_VEC_HPP

#include <cassert>
#include <complex>

namespace Utilities {

/** Very boring field element, could be used for transverse simulations. */
template <typename T>
class Vec1t
{
public:
  Vec1t(T _x=0.0) : x(_x) {}
  T x;

  /** Compute length of vector. */
  T length(void) const { return x; }

  /** Compute length^2 of vector. */
  T length2(void) const { return x*x; }

  /** Calculate the dot product. */
  T dot(const Vec1t<T>& a) { return x*a.x; }

  /** Number of components in vector. */
  constexpr static int veclength(void) { return 1; }

  Vec1t<T> normalized(void) const
  {
    return *this / length();
  }
};

/** 2-dimensional field element. */
template <typename T>
class Vec2t
{
public:
  Vec2t(T _x=0.0, T _y=0.0) : x(_x), y(_y) {}
  T x, y;

  /** Compute length of vector. */
  T length(void) const { return pow(x*x + y*y, 0.5); }

  /** Compute length^2 of vector. */
  T length2(void) const { return x*x + y*y; }

  /** Calculate the dot product. */
  T dot(const Vec2t<T>& a) { return x*a.x + y*a.y; }

  /** Number of components in vector. */
  constexpr static int veclength(void) { return 2; }

  Vec2t<T> normalized(void) const
  {
    return *this / length();
  }
};

/** 3-dimensional field element, the most common type. */
template <typename T>
class Vec3t
{
public:
  constexpr Vec3t(T _x=0.0, T _y=0.0, T _z=0.0)
      : x(_x), y(_y), z(_z) {}
  T x, y, z;

  template <typename O>
  Vec3t(const Vec3t<O>& o)
      : x(o.x), y(o.y), z(o.z) {}

  /** Negate the vector. */
  Vec3t<T> operator-(void)
  {
    return Vec3t<T>(-x, -y, -z);
  }

  /** Compute length of vector. */
  T length(void) const { return pow(x*x + y*y + z*z, 0.5); }

  /** Compute length^2 of vector. */
  T length2(void) const { return x*x + y*y + z*z; }

  /** Compute cross product of this and o. */
  Vec3t<T> cross(const Vec3t<T>& o) const
  {
    return Vec3t<T>(y*o.z - z*o.y, z*o.x - x*o.z, x*o.y - y*o.x);
  }

  /** Calculate the dot product. */
  T dot(const Vec3t<T>& a) { return x*a.x + y*a.y + z*a.z; }

  Vec3t<T> normalized(void) const
  {
    return *this / length();
  }

  T& operator[] (unsigned i) { return axis(i); }
  const T& operator[] (unsigned i) const { return axis(i); }

  /** Return a rotated copy of this vector. */
  Vec3t<T> rotated(double t, int axis) const
  {
    if (axis == 2)
    {
      return Vec3t<T>(x, y*cos(t) - z*sin(t), y*sin(t) + z*cos(t));
    }

    if (axis == 1)
    {
      return Vec3t<T>(x*cos(t) + z*sin(t), y, z*cos(t) - x*sin(t));
    }

    if (axis == 0)
    {
      return Vec3t<T>(x*cos(t) - y*sin(t), x*sin(t) + y*cos(t), z);
    }

    throw std::logic_error("Vec: Invalid axis for rotation.");
  }

  /** Retrieve component a. */
  constexpr T& axis(unsigned a)
  {
    // TODO: This feels like kludge
    if (a == 2) return x;
    if (a == 1) return y;
    if (a == 0) return z;
    assert(false);
  }

  /** Retrieve component a. */
  constexpr const T& axis(unsigned a) const
  {
    // TODO: This feels like kludge
    if (a == 2) return x;
    if (a == 1) return y;
    if (a == 0) return z;
    assert(false);
  }

  /** Retrieve a not component. */
  constexpr const T& naxis(int a, int b=-1) const
  {
    // TODO: This feels like kludge
    if (a != 2 && b != 2) return x;
    if (a != 1 && b != 1) return y;
    if (a != 0 && b != 0) return z;
    assert(false);
  }

  /** Retrieve a not component. */
  constexpr T& naxis(int a, int b=-1)
  {
    // TODO: This feels like kludge
    if (a != 2 && b != 2) return x;
    if (a != 1 && b != 1) return y;
    if (a != 0 && b != 0) return z;
    assert(false);
  }

  /** Calculate the absolute value of the vector. */
  Vec3t<decltype(std::abs(T(0.0)))> abs(void) const
  {
    return Vec3t<decltype(std::abs(T(0.0)))>
        (std::abs(x), std::abs(y), std::abs(z));
  }

  /** Calculate the absolute value of the vector. */
  Vec3t<decltype(std::real(T(0.0)))> real(void) const
  {
    return Vec3t<decltype(std::real(T(0.0)))>
        (std::real(x), std::real(y), std::real(z));
  }

  /** Calculate the absolute value of the vector. */
  Vec3t<decltype(std::imag(T(0.0)))> imag(void) const
  {
    return Vec3t<decltype(std::imag(T(0.0)))>
        (std::imag(x), std::imag(y), std::imag(z));
  }

  /** Calculate the complex conjugate of the vector. */
  Vec3t<T> conj(void) const
  {
    return Vec3t<T>(std::conj(x), std::conj(y), std::conj(z));
  }

  /** Number of components in vector. */
  constexpr static int veclength(void) { return 3; }
};

template <class O, class T>
Vec3t<T> operator* (const O& a, const Vec3t<T>& b)
{
  return Vec3t<T>(a*b.x, a*b.y, a*b.z);
}

template <class O, class T>
Vec3t<T> operator/ (const O& a, const Vec3t<T>& b)
{
  return Vec3t<T>(a/b.x, a/b.y, a/b.z);
}

template <class O, class T>
Vec3t<T> operator* (const Vec3t<T>& b, const O& a)
{
  return a * b;
}

template <class O, class T>
Vec3t<T> operator/ (const Vec3t<T>& b, const O& a)
{
  return Vec3t<T>(b.x/a, b.y/a, b.z/a);
}

template <class T>
Vec3t<T> operator- (const Vec3t<T>& a, const Vec3t<T>& b)
{
  return Vec3t<T>(a.x-b.x, a.y-b.y, a.z-b.z);
}

template <class T>
Vec3t<T> operator- (const Vec3t<T>& a, double b)
{
  return Vec3t<T>(a.x - b, a.y - b, a.z - b);
}

template <class T>
Vec3t<T> operator+ (const Vec3t<T>& a, const Vec3t<T>& b)
{
  return Vec3t<T>(a.x+b.x, a.y+b.y, a.z+b.z);
}

template <class T>
Vec3t<T> operator* (const Vec3t<T>& a, const Vec3t<T>& b)
{
  return Vec3t<T>(a.x*b.x, a.y*b.y, a.z*b.z);
}

template <class T>
Vec3t<T>& operator+= (Vec3t<T>& a, const Vec3t<T>& b)
{
  a.x += b.x; a.y += b.y; a.z += b.z;
  return a;
}

template <class T>
Vec3t<T>& operator-= (Vec3t<T>& a, const Vec3t<T>& b)
{
  a.x -= b.x; a.y -= b.y; a.z -= b.z;
  return a;
}

template <class T>
Vec3t<T>& operator*= (Vec3t<T>& a, const Vec3t<T>& b)
{
  a.x *= b.x; a.y *= b.y; a.z *= b.z;
  return a;
}

template <class T, class O>
Vec3t<T>& operator*= (Vec3t<T>& a, const O& b)
{
  a.x *= b; a.y *= b; a.z *= b;
  return a;
}

template <class T, class O>
Vec3t<T>& operator/= (Vec3t<T>& a, const O& b)
{
  a.x /= b; a.y /= b; a.z /= b;
  return a;
}

template <class T>
Vec3t<T> real(const Vec3t< std::complex<T> >& a)
{
  return Vec3t<T>(std::real(a.x), std::real(a.y), std::real(a.z));
}

template <class O, class T>
O& operator<<(O& ostream, const Vec3t<T>& a)
{
  ostream << a.x << '\t' << a.y << '\t' << a.z;
  return ostream;
}

typedef Vec1t<double> Vec1d;
typedef Vec1t<std::complex<double> > Vec1c;

typedef Vec2t<double> Vec2d;
typedef Vec2t<std::complex<double> > Vec2c;

typedef Vec3t<double> Vec3d;
typedef Vec3t<std::complex<double> > Vec3c;

} /* namespace Fdtd */

#endif /* #ifndef FDTD_VEC_HPP */

