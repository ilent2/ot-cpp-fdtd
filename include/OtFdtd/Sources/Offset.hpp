/* Sources/Offset.hpp - Coordinate offsets for beam sources.
 *
 * Copyright (C) 2016 Isaac Lenton (aka ilent2)
 */

#include "Utilities/Vec.hpp"

#ifndef FDTD_SOURCES_OFFSET_HPP
#define FDTD_SOURCES_OFFSET_HPP

namespace Fdtd {
namespace Offset {

class _ModifyDummy {};

/** Helpers for modifying coordinates. */
template <typename First=_ModifyDummy, typename... Modifiers>
struct Modify
{
  template <typename... Base>
  static Utilities::Vec3d modify(Utilities::Vec3d it, const Base&... b)
  {
    return Modify<Modifiers...>::modify(First::modify(it, b...), b...);
  }
};

/** Helpers for modifying coordinates. */
template <>
struct Modify<_ModifyDummy>
{
  template <typename... Base>
  static Utilities::Vec3d modify(Utilities::Vec3d it, const Base&... b)
  {
    return it;
  }
};

/** Subtract the coordinate system centre from the coordinate. */
class Centre
{
public:
  template <typename Base>
  static Utilities::Vec3d modify(Utilities::Vec3d it, const Base& base)
  {
    return it - base.centre();
  }

  static Utilities::Vec3d modify(Utilities::Vec3d it)
  {
    return Utilities::Vec3d(0.0, 0.0, 0.0);
  }
};

template <const Utilities::Vec3d& offset>
class Location
{
public:
  template <typename... Base>
  static Utilities::Vec3d modify(Utilities::Vec3d it, const Base&... b)
  {
    return it - offset;
  }
};

/** Subtract the location of the specified element.
 *
 * @tparam loc      The indices of the element.
 */
template <unsigned... loc>
class ElementLocation
{
public:
  template <typename T, typename Base>
  static T modify(const T& it, const Base& b)
  {
    return it - b.location(b.at(loc...));
  }
};

/** Apply a rotation to the coordinate around axis.
 *
 * @tparam Angle    The angle in radians.
 * @tparam Axis     The axis to rotate the vector around.
 */
template <const double& Angle, int Axis>
class Rotation
{
public:
  template <typename T, typename... Base>
  static T modify(const T& it, const Base&... b)
  {
    return it.rotated(Angle, Axis);
  }
};

} // namespace Offset
} // namespace Fdtd

#endif // #ifndef FDTD_SOURCES_OFFSET_HPP

