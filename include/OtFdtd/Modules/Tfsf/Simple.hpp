/* Fdtd/Tfsf/Simple.hpp - Automatically deduced TFSF.
 *
 * TODO: Generalize for 1-D and 2-D again.
 *
 * Copyright (C) 2016 Isaac Lenton (aka ilent2)
 */

#ifndef FDTD_TFSF_SIMPLE_HPP
#define FDTD_TFSF_SIMPLE_HPP

#include <vector>
#include "../EmptyModule.hpp"
#include "../../Utilities/Vec.hpp"

namespace Fdtd {
namespace Tfsf {

template <typename Source, unsigned Depth>
struct Simple : public EmptyModule
{

  /** Hook into the assembled object for returning the source power. */
  template <typename Base>
  static double sourcePower(Base& base)
  {
    // It doesn't make sense to have two TFSF sources of the same
    // type in the same location, so the following is safe.
    return base.template tfsfSourcePower<Source, Depth>();
  }

  template <typename Base> class EvolveMethods : public Base
  {
    Source m_source;
    double m_sourcePower;
    bool m_bRecalculatePower;

    static constexpr unsigned imin = Depth;
    static constexpr unsigned imax = Base::length(2) - Depth - 1;
    static constexpr unsigned jmin = Depth;
    static constexpr unsigned jmax = Base::length(1) - Depth - 1;
    static constexpr unsigned kmin = Depth;
    static constexpr unsigned kmax = Base::length(0) - Depth - 1;

    /** Create a vector of locations the E-field needs. */
    std::vector<Utilities::Vec3d> hfieldLocations(void)
    {
      std::vector<Utilities::Vec3d> ret;

      for (unsigned j = jmin+1; j < jmax; ++j)
      {
        for (unsigned k = kmin+1; k < kmax; ++k)
        {
          ret.push_back(Base::hcoord(Base::at(imin+1, j, k), 0));
          ret.push_back(Base::hcoord(Base::at(imin+1, j, k), 1));
          ret.push_back(Base::hcoord(Base::at(imax, j, k), 0));
          ret.push_back(Base::hcoord(Base::at(imax, j, k), 1));
        }
      }

      for (unsigned i = imin+1; i < imax; ++i)
      {
        for (unsigned j = jmin+1; j < jmax; ++j)
        {
          ret.push_back(Base::hcoord(Base::at(i, j, kmin+1), 1));
          ret.push_back(Base::hcoord(Base::at(i, j, kmin+1), 2));
          ret.push_back(Base::hcoord(Base::at(i, j, kmax), 1));
          ret.push_back(Base::hcoord(Base::at(i, j, kmax), 2));
        }

        for (unsigned k = kmin+1; k < kmax; ++k)
        {
          ret.push_back(Base::hcoord(Base::at(i, jmin+1, k), 0));
          ret.push_back(Base::hcoord(Base::at(i, jmin+1, k), 2));
          ret.push_back(Base::hcoord(Base::at(i, jmax, k), 0));
          ret.push_back(Base::hcoord(Base::at(i, jmax, k), 2));
        }
      }

      return ret;
    }

    /** Create a vector of locations the H-field needs. */
    std::vector<Utilities::Vec3d> efieldLocations(void)
    {
      std::vector<Utilities::Vec3d> ret;

      for (unsigned j = jmin+1; j < jmax; ++j)
      {
        for (unsigned k = kmin+1; k < kmax; ++k)
        {
          ret.push_back(Base::ecoord(Base::at(imin, j, k), 0));
          ret.push_back(Base::ecoord(Base::at(imin, j, k), 1));
          ret.push_back(Base::ecoord(Base::at(imax-1, j, k), 0));
          ret.push_back(Base::ecoord(Base::at(imax-1, j, k), 1));
        }
      }

      for (unsigned i = imin+1; i < imax; ++i)
      {
        for (unsigned j = jmin+1; j < jmax; ++j)
        {
          ret.push_back(Base::ecoord(Base::at(i, j, kmin), 1));
          ret.push_back(Base::ecoord(Base::at(i, j, kmin), 2));
          ret.push_back(Base::ecoord(Base::at(i, j, kmax-1), 1));
          ret.push_back(Base::ecoord(Base::at(i, j, kmax-1), 2));
        }

        for (unsigned k = kmin+1; k < kmax; ++k)
        {
          ret.push_back(Base::ecoord(Base::at(i, jmin, k), 0));
          ret.push_back(Base::ecoord(Base::at(i, jmin, k), 2));
          ret.push_back(Base::ecoord(Base::at(i, jmax-1, k), 0));
          ret.push_back(Base::ecoord(Base::at(i, jmax-1, k), 2));
        }
      }

      return ret;
    }

    /** Create a vector of locations for the power calculation. */
    std::vector<Utilities::Vec3d> powerLocations(void)
    {
      std::vector<Utilities::Vec3d> ret;

      for (unsigned j = jmin+1; j < jmax; ++j)
      {
        for (unsigned k = kmin+1; k < kmax; ++k)
        {
          ret.push_back(Base::ecoord(Base::at(imin, j, k), -1));
          ret.push_back(Base::ecoord(Base::at(imax, j, k), -1));
        }
      }

      for (unsigned i = imin+1; i < imax; ++i)
      {
        for (unsigned j = jmin+1; j < jmax; ++j)
        {
          ret.push_back(Base::ecoord(Base::at(i, j, kmin), -1));
          ret.push_back(Base::ecoord(Base::at(i, j, kmax), -1));
        }

        for (unsigned k = kmin+1; k < kmax; ++k)
        {
          ret.push_back(Base::ecoord(Base::at(i, jmin, k), -1));
          ret.push_back(Base::ecoord(Base::at(i, jmax, k), -1));
        }
      }

      return ret;
    }

  public:

    /** Report information about us and our source. */
    template <typename T> void report_initialize(T& fp) const
    {
      fp << "Using Tfsf::Simple at depth " << Depth << " with..." << std::endl;
      m_source.report_initialize(fp);
      Base::report_initialize(fp);
    }

    /** Initialize the source. */
    void initialize(void)
    {
      Base::initialize();
      m_source.initialize(*this, efieldLocations(), hfieldLocations());

      // Calculate the source power
      m_bRecalculatePower = true;
    }

    /** Return the source power. */
    template <typename S, unsigned D> std::enable_if_t<
        std::is_same<S, Source>::value && D == Depth, double>
    tfsfSourcePower(void)
    {
      if (m_bRecalculatePower)
      {
        m_sourcePower = calculateSourcePower();
        m_bRecalculatePower = false;
      }

      return m_sourcePower;
    }

    /** Calculate the source power. */
    double calculateSourcePower(void)
    {
      double power = 0.0;

      Source powerSource;
      powerSource.initialize(*this, powerLocations());
      auto vecEData = powerSource.cefield();
      auto vecHData = powerSource.chfield();

      auto hData = vecHData.cbegin();
      auto eData = vecEData.cbegin();

      for (unsigned j = jmin+1; j < jmax; ++j)
      {
        for (unsigned k = kmin+1; k < kmax; ++k)
        {
          auto em0 = Base::at(imin, j, k);
          auto em1 = Base::at(imax, j, k);

          Utilities::Vec3d norm0, norm1;
          norm0.axis(2) = Base::earea(em0, 2);
          norm1.axis(2) = Base::earea(em1, 2);

          double imp0 = pow(em0->permeability()/em0->isoPermittivity(), 0.5);
          double imp1 = pow(em1->permeability()/em1->isoPermittivity(), 0.5);

          auto h0 = (*hData++).conj()/imp0;
          auto h1 = (*hData++).conj()/imp1;
          auto e0 = *eData++;
          auto e1 = *eData++;

          power += norm0.dot(e0.cross(h0).abs());
          power += norm1.dot(e1.cross(h1).abs());
        }
      }

      for (unsigned i = imin+1; i < imax; ++i)
      {
        for (unsigned j = jmin+1; j < jmax; ++j)
        {
          auto em0 = Base::at(i, j, kmin);
          auto em1 = Base::at(i, j, kmax);

          Utilities::Vec3d norm0, norm1;
          norm0.axis(0) = Base::earea(em0, 0);
          norm1.axis(0) = Base::earea(em1, 0);

          double imp0 = pow(em0->permeability()/em0->isoPermittivity(), 0.5);
          double imp1 = pow(em1->permeability()/em1->isoPermittivity(), 0.5);

          auto h0 = (*hData++).conj()/imp0;
          auto h1 = (*hData++).conj()/imp1;
          auto e0 = *eData++;
          auto e1 = *eData++;

          power += norm0.dot(e0.cross(h0).abs());
          power += norm1.dot(e1.cross(h1).abs());
        }

        for (unsigned k = kmin+1; k < kmax; ++k)
        {
          auto em0 = Base::at(i, jmin, k);
          auto em1 = Base::at(i, jmax, k);

          Utilities::Vec3d norm0, norm1;
          norm0.axis(1) = Base::earea(em0, 1);
          norm1.axis(1) = Base::earea(em1, 1);

          double imp0 = pow(em0->permeability()/em0->isoPermittivity(), 0.5);
          double imp1 = pow(em1->permeability()/em1->isoPermittivity(), 0.5);

          auto h0 = (*hData++).conj()/imp0;
          auto h1 = (*hData++).conj()/imp1;
          auto e0 = *eData++;
          auto e1 = *eData++;

          power += norm0.dot(e0.cross(h0).abs());
          power += norm1.dot(e1.cross(h1).abs());
        }
      }

      return power * 0.5;
    }

  // TODO: Another approach could be to store the element update
  //    equation with the element, rather than have multiple
  //    advanceE/advanceH methods.
  //
  //    If we are processor limited this method might be faster
  //    but more memory intensive.

  // TODO: An another approach: Have a single advanceE method,
  //    but implement AdvanceMethods:: member functions that
  //    determine if they should act on a particular element.
  //
  //    This method would rely on the compiler to make optimisations,
  //    the net result could be the same as the above method.

    /** The 1-step update equation for the E-field. */
    void advanceE(typename Base::TimingIterator it)
    {
      auto vecTfsfData = m_source.hfield(it->hfieldT());
      auto tfsfData = vecTfsfData.cbegin();

      for (unsigned j = jmin+1; j < jmax; ++j)
      {
        for (unsigned k = kmin+1; k < kmax; ++k)
        {
          // Update terms with x-derivative
          auto em0 = Base::at(imin, j, k);
          auto em1 = Base::at(imax-1, j, k);

          auto h0z = *tfsfData++;
          auto h0y = *tfsfData++;
          auto h1z = *tfsfData++;
          auto h1y = *tfsfData++;

          // TODO: An impedance member would be useful
          //double eps1 = Base::hstep(em1, 2)*Base::impedance(em1);

          double dx0 = Base::hstep(em0, 2);
          double dx1 = Base::hstep(em1, 2);
          double eps0 = dx0*pow(em0->permeability()/em0->isoPermittivity(), 0.5);
          double eps1 = dx1*pow(em1->permeability()/em1->isoPermittivity(), 0.5);

          Base::dfieldAdd(em0, -it->efieldDt()/eps0*h0z.z, 1);
          Base::dfieldAdd(em0,  it->efieldDt()/eps0*h0y.y, 0);
          Base::dfieldAdd(em1,  it->efieldDt()/eps1*h1z.z, 1);
          Base::dfieldAdd(em1, -it->efieldDt()/eps1*h1y.y, 0);
        }
      }

      for (unsigned i = imin+1; i < imax; ++i)
      {
        for (unsigned j = jmin+1; j < jmax; ++j)
        {
          // Update terms with z-derivative
          auto em0 = Base::at(i, j, kmin);
          auto em1 = Base::at(i, j, kmax-1);

          auto h0y = *tfsfData++;
          auto h0x = *tfsfData++;
          auto h1y = *tfsfData++;
          auto h1x = *tfsfData++;

          double dz0 = Base::hstep(em0, 0);
          double dz1 = Base::hstep(em1, 0);
          double eps0 = dz0*pow(em0->permeability()/em0->isoPermittivity(), 0.5);
          double eps1 = dz1*pow(em1->permeability()/em1->isoPermittivity(), 0.5);

          Base::dfieldAdd(em0, -it->efieldDt()/eps0*h0y.y, 2);
          Base::dfieldAdd(em0,  it->efieldDt()/eps0*h0x.x, 1);
          Base::dfieldAdd(em1,  it->efieldDt()/eps1*h1y.y, 2);
          Base::dfieldAdd(em1, -it->efieldDt()/eps1*h1x.x, 1);
        }

        for (unsigned k = kmin+1; k < kmax; ++k)
        {
          // Update terms with y-derivative
          auto em0 = Base::at(i, jmin, k);
          auto em1 = Base::at(i, jmax-1, k);

          auto h0z = *tfsfData++;
          auto h0x = *tfsfData++;
          auto h1z = *tfsfData++;
          auto h1x = *tfsfData++;

          double dy0 = Base::hstep(em0, 1);
          double dy1 = Base::hstep(em1, 1);
          double eps0 = dy0*pow(em0->permeability()/em0->isoPermittivity(), 0.5);
          double eps1 = dy1*pow(em1->permeability()/em1->isoPermittivity(), 0.5);

          Base::dfieldAdd(em0,  it->efieldDt()/eps0*h0z.z, 2);
          Base::dfieldAdd(em0, -it->efieldDt()/eps0*h0x.x, 0);
          Base::dfieldAdd(em1, -it->efieldDt()/eps1*h1z.z, 2);
          Base::dfieldAdd(em1,  it->efieldDt()/eps1*h1x.x, 0);
        }
      }

      Base::advanceE(it);
    }

    /** The 1-step update equation for the H-field. */
    void advanceH(typename Base::TimingIterator it)
    {
      auto vecTfsfData = m_source.efield(it->efieldTlast());
      auto tfsfData = vecTfsfData.cbegin();

      for (unsigned j = jmin+1; j < jmax; ++j)
      {
        for (unsigned k = kmin+1; k < kmax; ++k)
        {
          // Update terms with x-derivative
          auto em0 = Base::at(imin+1, j, k);
          auto em1 = Base::at(imax, j, k);

          auto e0z = *tfsfData++;
          auto e0y = *tfsfData++;
          auto e1z = *tfsfData++;
          auto e1y = *tfsfData++;

          double dx0 = Base::estep(em0, 2);
          double dx1 = Base::estep(em1, 2);

          Base::bfieldAdd(em0, -it->hfieldDt()/dx0*e0z.z, 1);
          Base::bfieldAdd(em0,  it->hfieldDt()/dx0*e0y.y, 0);
          Base::bfieldAdd(em1,  it->hfieldDt()/dx1*e1z.z, 1);
          Base::bfieldAdd(em1, -it->hfieldDt()/dx1*e1y.y, 0);
        }
      }

      for (unsigned i = imin+1; i < imax; ++i)
      {
        for (unsigned j = jmin+1; j < jmax; ++j)
        {
          // Update terms with z-derivative
          auto em0 = Base::at(i, j, kmin+1);
          auto em1 = Base::at(i, j, kmax);

          auto e0y = *tfsfData++;
          auto e0x = *tfsfData++;
          auto e1y = *tfsfData++;
          auto e1x = *tfsfData++;

          double dz0 = Base::estep(em0, 0);
          double dz1 = Base::estep(em1, 0);

          Base::bfieldAdd(em0, -it->hfieldDt()/dz0*e0y.y, 2);
          Base::bfieldAdd(em0,  it->hfieldDt()/dz0*e0x.x, 1);
          Base::bfieldAdd(em1,  it->hfieldDt()/dz1*e1y.y, 2);
          Base::bfieldAdd(em1, -it->hfieldDt()/dz1*e1x.x, 1);
        }

        for (unsigned k = kmin+1; k < kmax; ++k)
        {
          // Update terms with y-derivative
          auto em0 = Base::at(i, jmin+1, k);
          auto em1 = Base::at(i, jmax, k);

          auto e0z = *tfsfData++;
          auto e0x = *tfsfData++;
          auto e1z = *tfsfData++;
          auto e1x = *tfsfData++;

          double dy0 = Base::estep(em0, 1);
          double dy1 = Base::estep(em1, 1);

          Base::bfieldAdd(em0,  it->hfieldDt()/dy0*e0z.z, 2);
          Base::bfieldAdd(em0, -it->hfieldDt()/dy0*e0x.x, 0);
          Base::bfieldAdd(em1, -it->hfieldDt()/dy1*e1z.z, 2);
          Base::bfieldAdd(em1,  it->hfieldDt()/dy1*e1x.x, 0);
        }
      }

      Base::advanceH(it);
    }
  };

};

} // namespace Tfsf
} // namespace Fdtd

#endif // #ifndef FDTD_TFSF_SIMPLE_HPP

