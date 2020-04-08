/* Coordinate/Stretched.hpp - A stretched coordinate type.
 *
 * This type requires an existing axis to stretch.
 *
 * Copyright (C) 2016 Isaac Lenton (aka ilent2)
 */

#ifndef FDTD_COORDINATE_STRETCHED_HPP
#define FDTD_COORDINATE_STRETCHED_HPP

#include <array>
#include "../EmptyModule.hpp"

namespace Fdtd {
namespace Coordinate {

/** A stretched coordinate axis.
 *
 * @tparam Index      : The index of the axis to stretch.
 * @tparam Grading    : A type with estep and hstep methods.
 *
 *    The e/hstep methods should take three arguments, a reference
 *    to the Simulation data, a linear index and a axis.
 */
template <int Index, typename Grading>
struct Stretched : public EmptyModule
{

  template <typename Base> class CoordinateMethods : public Base
  {
    // Check that Index is valid
    static_assert(Index < Base::CoordinateIndex,
        "Invalid coordinate axis selected for stretching.");

    std::array<double, Base::length(Index)> m_estep;
    std::array<double, Base::length(Index)> m_hstep;
    std::array<double, Base::length(Index)> m_elength;
    std::array<double, Base::length(Index)> m_hlength;

  public:

    /** Report we are using linear coordinates on this axis. */
    template <typename T> void report_construct(T& fp) const
    {
      fp << "Applying Coordinates::Stretched to axis " << Index << std::endl;
      Base::report_construct(fp);
    }

    /** Pre-compute the locations and step sizes. */
    void initialize(void)
    {
      Base::initialize();

      auto loc = Base::dataBegin();

      m_elength[0] = Grading::offset(*this, Index, Base::elength(loc, Index));
      m_hlength[0] = 0.0;

      for (unsigned i = 0; i < Base::length(Index); ++i)
      {
        m_estep[i] = Base::estep(loc, Index) * Grading::estep(*this, i, Index);
        m_hstep[i] = Base::hstep(loc, Index) * Grading::hstep(*this, i, Index);

        if (i > 0)
        {
          m_elength[i] = m_elength[i-1] + m_estep[i];
          m_hlength[i] = m_hlength[i-1] + m_hstep[i-1];
        }

        loc = Base::next(loc, Index);
      }
    }

    /** Return the distance between two adjacent H-field locations. */
    double hstep(typename Base::DataIterator it,
        int Axis=Index) const
    {
      if (Axis != Index) return Base::hstep(it, Axis);
      return m_hstep[Base::getOffset(it, Axis)];
    }

    /** Return the distance between two adjacent E-field locations. */
    double estep(typename Base::DataIterator it,
        int Axis=Index) const
    {
      if (Axis != Index) return Base::estep(it, Axis);
      return m_estep[Base::getOffset(it, Axis)];
    }

    /** Compute the offset of the component vector from the origin. */
    double elength(typename Base::DataIterator it,
        int Axis=Index, int Component=-1) const
    {
      if (Axis != Index) return Base::elength(it, Axis, Component);
      unsigned i = Base::getOffset(it, Axis);
      return Component == Axis ? m_hlength[i] : m_elength[i];
    }

    /** Compute the offset of the component vector from the origin. */
    double hlength(typename Base::DataIterator it,
        int Axis=Index, int Component=-1) const
    {
      if (Axis != Index) return Base::hlength(it, Axis, Component);
      unsigned i = Base::getOffset(it, Axis);
      return Component == Axis ? m_elength[i] : m_hlength[i];
    }
  };
  
};

} // namespace Coordinate
} // namespace Fdtd

#endif // #ifndef FDTD_COORDINATE_STRETCHED_HPP

