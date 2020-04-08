/* Output/Integrate/Surface.hpp - A class describing the surface integral.
 *
 * Copyright (C) 2016 Isaac Lenton (aka ilent2)
 */

#ifndef FDTD_OUTPUT_INTEGRATE_ITERATORS_HPP
#define FDTD_OUTPUT_INTEGRATE_ITERATORS_HPP

namespace Fdtd {
namespace Output {
namespace Integrate {
namespace Iterators {

/** Provides the iteration method for a single plane.
 *
 * @tparam Depth    : Description of the output surface.
 * @tparam Axis     : The axis normal to the integration plane.
 */
template <typename Depth, int Axis>
struct SurfacePlane
{
  static std::string repr(void)
  {
    return "SurfacePlain on axis " + std::to_string(Axis)
        + " with depth " + Depth::repr();
  }

  template <typename B, typename Kernel>
  struct Iterator : public Kernel::template Parameter<B>
  {
    using Base = typename Kernel::template Parameter<B>;

    static_assert(Axis >= 0 && Axis < 3, "Only implemented for up to 3-D");

    static_assert(Depth::x0 + Depth::x1 < Base::length(2),
        "Total Depth::x[01] should be less than length(2).");
    static_assert(Depth::y0 + Depth::y1 < Base::length(1),
        "Total Depth::y[01] should be less than length(1).");
    static_assert(Depth::z0 + Depth::z1 < Base::length(0),
        "Total Depth::z[01] should be less than length(0).");

  public:

    Utilities::Vec3d Integrate(typename Base::TimingIterator it)
    {
      Utilities::Vec3d ret(0.0, 0.0, 0.0);

      unsigned imin = Depth::x0;
      unsigned imax = Base::length(2) - Depth::x1;
      unsigned jmin = Depth::y0;
      unsigned jmax = Base::length(1) - Depth::y1;
      unsigned kmin = Depth::z0;
      unsigned kmax = Base::length(0) - Depth::z1;

      if (Axis == 0)
      {
        for (unsigned i = imin; i < imax; ++i)
        {
          for (unsigned j = jmin; j < jmax; ++j)
          {
            auto loc = Base::at(i, j, (kmin + kmax)/2);
            Utilities::Vec3d norm; norm.axis(Axis) = Base::earea(loc, Axis);
            ret += Base::SurfaceKernel(it, loc, norm, Axis);
          }
        }
      }
      else if (Axis == 1)
      {
        for (unsigned i = imin; i < imax; ++i)
        {
          for (unsigned k = kmin; k < kmax; ++k)
          {
            auto loc = Base::at(i, (jmin + jmax)/2, k);
            Utilities::Vec3d norm; norm.axis(Axis) = Base::earea(loc, Axis);
            ret += Base::SurfaceKernel(it, loc, norm, Axis);
          }
        }
      }
      else if (Axis == 2)
      {
        for (unsigned j = jmin; j < jmax; ++j)
        {
          for (unsigned k = kmin; k < kmax; ++k)
          {
            auto loc = Base::at((imin + imax)/2, j, k);
            Utilities::Vec3d norm; norm.axis(Axis) = Base::earea(loc, Axis);
            ret += Base::SurfaceKernel(it, loc, norm, Axis);
          }
        }
      }

      return ret;
    }
  };
};

/** Provides the iteration method for a pair of parallel planes.
 *
 * @tparam Base     : Base class to be decorated.
 * @tparam Depth    : Description of the output surface.
 * @tparam Axis     : The axis normal to the integration plane.
 * @tparam Kernel   : The kernel to use for each integration element.
 *
 * TODO: This might be better with template specialisations of Axis
 */
template <typename Depth, int Axis>
struct SurfacePair
{
  static std::string repr(void)
  {
    return "SurfacePair on axis " + std::to_string(Axis)
        + " with depth " + Depth::repr();
  }

  template <typename B, typename Kernel>
  class Iterator : public Kernel::template Parameter<B>
  {
    using Base = typename Kernel::template Parameter<B>;

    static_assert(Axis >= 0 && Axis < 3, "Only implemented for up to 3-D");

    static_assert(Depth::x0 + Depth::x1 < Base::length(2),
        "Total Depth::x[01] should be less than length(2).");
    static_assert(Depth::y0 + Depth::y1 < Base::length(1),
        "Total Depth::y[01] should be less than length(1).");
    static_assert(Depth::z0 + Depth::z1 < Base::length(0),
        "Total Depth::z[01] should be less than length(0).");

  public:

    Utilities::Vec3d Integrate(typename Base::TimingIterator it)
    {
      Utilities::Vec3d ret(0.0, 0.0, 0.0);

      unsigned imin = Depth::x0;
      unsigned imax = Base::length(2) - Depth::x1;
      unsigned jmin = Depth::y0;
      unsigned jmax = Base::length(1) - Depth::y1;
      unsigned kmin = Depth::z0;
      unsigned kmax = Base::length(0) - Depth::z1;

      if (Axis == 0)
      {
        for (unsigned i = imin; i < imax; ++i)
        {
          for (unsigned j = jmin; j < jmax; ++j)
          {
            auto loct = Base::at(i, j, kmin);
            Utilities::Vec3d normt; normt.axis(Axis) = Base::earea(loct, Axis);

            auto locb = Base::at(i, j, kmax);
            Utilities::Vec3d normb; normb.axis(Axis) = Base::earea(locb, Axis);

            ret += Base::SurfaceKernel(it, loct, normt, Axis);
            ret += Base::SurfaceKernel(it, locb, -normb, Axis);
          }
        }
      }
      else if (Axis == 1)
      {
        // TODO: Axis == 1 || 2 could probably be evaluate more efficiently.
        for (unsigned i = imin; i < imax; ++i)
        {
          for (unsigned k = kmin; k < kmax; ++k)
          {
            auto loct = Base::at(i, jmin, k);
            Utilities::Vec3d normt; normt.axis(Axis) = Base::earea(loct, Axis);

            auto locb = Base::at(i, jmax, k);
            Utilities::Vec3d normb; normb.axis(Axis) = Base::earea(locb, Axis);

            ret += Base::SurfaceKernel(it, loct, normt, Axis);
            ret += Base::SurfaceKernel(it, locb, -normb, Axis);
          }
        }
      }
      else if (Axis == 2)
      {
        // TODO: Axis == 1 || 2 could probably be evaluate more efficiently.
        for (unsigned j = jmin; j < jmax; ++j)
        {
          for (unsigned k = kmin; k < kmax; ++k)
          {
            auto loct = Base::at(imin, j, k);
            Utilities::Vec3d normt; normt.axis(Axis) = Base::earea(loct, Axis);

            auto locb = Base::at(imax, j, k);
            Utilities::Vec3d normb; normb.axis(Axis) = Base::earea(locb, Axis);

            ret += Base::SurfaceKernel(it, loct, normt, Axis);
            ret += Base::SurfaceKernel(it, locb, -normb, Axis);
          }
        }
      }

      return -ret;
    }
  };
};

/** Provides the iteration method for a cube surface.
 *
 * @tparam Depth    : Description of the output surface.
 */
template <typename Depth>
struct Surface
{
  static std::string repr(void)
  {
    return "Surface with depth " + Depth::repr();
  }

  template <typename B, typename Kernel>
  class Iterator : public Kernel::template Parameter<B>
  {
    using Base = typename Kernel::template Parameter<B>;

    static_assert(Depth::x0 + Depth::x1 < Base::length(2),
        "Total Depth::x[01] should be less than length(2).");
    static_assert(Depth::y0 + Depth::y1 < Base::length(1),
        "Total Depth::y[01] should be less than length(1).");
    static_assert(Depth::z0 + Depth::z1 < Base::length(0),
        "Total Depth::z[01] should be less than length(0).");

  public:

    Utilities::Vec3d Integrate(typename Base::TimingIterator it)
    {
      Utilities::Vec3d ret(0.0, 0.0, 0.0);

      unsigned imin = Depth::x0;
      unsigned imax = Base::length(2) - Depth::x1;
      unsigned jmin = Depth::y0;
      unsigned jmax = Base::length(1) - Depth::y1;
      unsigned kmin = Depth::z0;
      unsigned kmax = Base::length(0) - Depth::z1;

      for (unsigned i = imin; i < imax; ++i)
      {
        for (unsigned j = jmin; j < jmax; ++j)
        {
          int Axis = 0;

          auto loct = Base::at(i, j, kmin);
          Utilities::Vec3d normt; normt.axis(Axis) = Base::earea(loct, Axis);

          auto locb = Base::at(i, j, kmax);
          Utilities::Vec3d normb; normb.axis(Axis) = Base::earea(locb, Axis);

          ret += Base::SurfaceKernel(it, loct, normt, Axis);
          ret += Base::SurfaceKernel(it, locb, -normb, Axis);
        }
      }

      for (unsigned i = imin; i < imax; ++i)
      {
        for (unsigned k = kmin; k < kmax; ++k)
        {
          int Axis = 1;

          auto loct = Base::at(i, jmin, k);
          Utilities::Vec3d normt; normt.axis(Axis) = Base::earea(loct, Axis);

          auto locb = Base::at(i, jmax, k);
          Utilities::Vec3d normb; normb.axis(Axis) = Base::earea(locb, Axis);

          ret += Base::SurfaceKernel(it, loct, normt, Axis);
          ret += Base::SurfaceKernel(it, locb, -normb, Axis);
        }
      }

      for (unsigned j = jmin; j < jmax; ++j)
      {
        for (unsigned k = kmin; k < kmax; ++k)
        {
          int Axis = 2;

          auto loct = Base::at(imin, j, k);
          Utilities::Vec3d normt; normt.axis(Axis) = Base::earea(loct, Axis);

          auto locb = Base::at(imax, j, k);
          Utilities::Vec3d normb; normb.axis(Axis) = Base::earea(locb, Axis);

          ret += Base::SurfaceKernel(it, loct, normt, Axis);
          ret += Base::SurfaceKernel(it, locb, -normb, Axis);
        }
      }

      return -ret;
    }
  };
};

/** Provides the iteration method for a cube volume.
 *
 * @tparam Depth    : Description of the output surface.
 */
template <typename Depth>
struct Volume
{
  static std::string repr(void)
  {
    return "Volume with depth " + Depth::repr();
  }

  template <typename B, typename Kernel>
  class Iterator : public Kernel::template Parameter<B>
  {
    using Base = typename Kernel::template Parameter<B>;

  public:

    Utilities::Vec3d Integrate(typename Base::TimingIterator it)
    {
      Utilities::Vec3d ret;

      unsigned imin = Depth::x0;
      unsigned imax = Base::length(2) - Depth::x1;
      unsigned jmin = Depth::y0;
      unsigned jmax = Base::length(1) - Depth::y1;
      unsigned kmin = Depth::z0;
      unsigned kmax = Base::length(0) - Depth::z1;

      for (unsigned i = imin; i < imax; ++i)
      {
        for (unsigned j = jmin; j < jmax; ++j)
        {
          for (unsigned k = kmin; k < kmax; ++k)
          {
            auto loc = Base::at(i, j, k);
            double volume = Base::estep(loc, 0)
                          * Base::estep(loc, 1)
                          * Base::estep(loc, 2);
            ret += Base::VolumeKernel(it, loc, volume);
          }
        }
      }

      return ret;
    }
  };
};

} // namespace Iterators
} // namespace Integrate
} // namespace Output
} // namespace Fdtd

#endif // #ifndef FDTD_OUTPUT_INTEGRATE_ITERATORS_HPP

