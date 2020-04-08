/* Cpml.hpp - CPML implementation for half of one axis.
 *
 * Copyright (C) 2016 Isaac Lenton (aka ilent2)
 */

#ifndef FDTD_CPML_CPML_HPP
#define FDTD_CPML_CPML_HPP

#include "Utilities/Vec.hpp"
#include <vector>

namespace Fdtd {
namespace Cpml {

/** Recursive application of CPML to the different axes. */
template <int axis, int Axis>
class CpmlAxesFor
{
  constexpr static int naxis = ~(axis ^ Axis) & 0x3;

public:
  /** Apply the CPML step to the E-field. */
  template <typename Base>
  static void cpmlE(Base& base, typename Base::TimingIterator it,
      unsigned depth, typename Base::DataIterator loc, unsigned count,
      double hstep)
  {
    if (axis != Axis)
    {
      const auto& h0 = loc->hfield();
      const auto h1 = base.nextr(loc, Axis)->hfield();

      double diff = (h1.axis(axis) - h0.axis(axis))/hstep;
      double phi = base.m_ephi[count].axis(axis)
          = base.m_eBeta[depth]*base.m_ephi[count].axis(axis)
          + base.m_eAlpha[depth]*diff;

      if ((Axis % 2) ^ (axis % 2) ^ (Axis < axis))
      {
        base.dfieldAdd(loc, it->hfieldDt() * phi, naxis);
      }
      else
      {
        base.dfieldAdd(loc, -it->hfieldDt() * phi, naxis);
      }
    }

    CpmlAxesFor<axis-1, Axis>::cpmlE(
        base, it, depth, loc, count, hstep);
  }

  /** Apply the CPML step to the H-field. */
  template <typename Base>
  static void cpmlH(Base& base, typename Base::TimingIterator it,
      unsigned depth, typename Base::DataIterator loc, unsigned count,
      double estep)
  {
    if (axis != Axis)
    {
      const auto e0 = base.prevr(loc, Axis)->efield();
      const auto& e1 = loc->efield();

      double diff = (e1.axis(axis) - e0.axis(axis))/estep;
      double phi = base.m_hphi[count].axis(axis)
          = base.m_hBeta[depth]*base.m_hphi[count].axis(axis)
          + base.m_hAlpha[depth]*diff;

      if ((Axis % 2) ^ (axis % 2) ^ (Axis < axis))
      {
        base.bfieldAdd(loc, -it->hfieldDt() * phi, naxis);
      }
      else
      {
        base.bfieldAdd(loc, it->hfieldDt() * phi, naxis);
      }
    }

    CpmlAxesFor<axis-1, Axis>::cpmlH(
        base, it, depth, loc, count, estep);
  }
};

/** Terminating case for CpmlAxes. */
template <int Axis>
class CpmlAxesFor<-1, Axis>
{
public:
  template <typename Base>
  static void cpmlE(Base& base, typename Base::TimingIterator it,
      unsigned depth, typename Base::DataIterator loc, unsigned count,
      double hstep) {}
  template <typename Base>
  static void cpmlH(Base& base, typename Base::TimingIterator it,
      unsigned depth, typename Base::DataIterator loc, unsigned count,
      double estep) {}
};


/** Implementation of CPML for one end of a axis.
 *
 * @tparam Base         Base class used by Fdtd::Fdtd.
 * @tparam Depth        Depth of CPML layer from simulation boundary.
 * @tparam Axis         Axis to add CPML to.
 * @tparam Forward      End to apply CPML to (true: start, false: end).
 * @tparam eps0         Permittivity inside CPML.
 * @tparam mu0          Permeability inside CPML.
 * @tparam Grading      Grading used for CPML parameters.
 */
template <typename Base, unsigned Depth, int Axis, bool Forward,
    const double& TimeStep, const double& eps0, const double& mu0,
    typename Grading = DefaultGrading>
class Cpml : public Base
{
  // TODO: I've hard coded veclength to 3, not sure how we are
  //    going to allow variable field element data types at the moment.
  template <int a, int c> friend class CpmlAxesFor;
  using CpmlAxes = CpmlAxesFor<3-1, Axis>;

  static_assert(Depth < Base::length(Axis),
      "Depth must be in the range of valid indices for Axis.");

  static_assert(Base::AxisIndex == 3,
      "This version only supports D=3, see CpmlGeneral");

  // Note: For now we are using Vec3d for m_ephi and m_hphi,
  //    but we only need a Vec2d, however axis/naxis...
  //  Hmm, Perhaps we should include an assertion on EmElement type.

  std::vector<Utilities::Vec3d> m_ephi;      ///< Convolution data for E-field
  std::vector<Utilities::Vec3d> m_hphi;      ///< Convolution data for H-field

  std::vector<double> m_eAlpha;   ///< Graded E-field alpha
  std::vector<double> m_hAlpha;   ///< Graded H-field alpha

  std::vector<double> m_eBeta;    ///< Graded E-field beta
  std::vector<double> m_hBeta;    ///< Graded H-field beta

  static constexpr unsigned amin =
      Forward ? 0 : Base::length(Axis) - Depth;
  static constexpr unsigned imin = Axis == 2 ? amin : 0;
  static constexpr unsigned jmin = Axis == 1 ? amin : 0;
  static constexpr unsigned kmin = Axis == 0 ? amin : 0;
  static constexpr unsigned imax = Axis == 2 ? amin + Depth : Base::length(2);
  static constexpr unsigned jmax = Axis == 1 ? amin + Depth : Base::length(1);
  static constexpr unsigned kmax = Axis == 0 ? amin + Depth : Base::length(0);

public:

  /** Allocate memory for convolution. */
  Cpml(void)
  {
    // Calculate the required amount of memory for convolution data
    unsigned phiLength = Depth;
    for (int i = 0; i < Base::AxisIndex; ++i)
    {
      if (i == Axis) continue;
      phiLength *= std::max(1u, Base::length(i));
    }

    // Request convolution data memory
    m_ephi.resize(phiLength);
    m_hphi.resize(phiLength);

    m_eAlpha.resize(Depth);
    m_hAlpha.resize(Depth);
    m_eBeta.resize(Depth);
    m_hBeta.resize(Depth);
  }

  /** Report our existence. */
  template <typename T>
  void report_construct(T& fp) const
  {
    fp << "Using Cpml on axis " << Axis << " at depth " << Depth
        << " with forward " << std::boolalpha <<  Forward << std::endl;
    Base::report_construct(fp);
  }

  /** Initialize memory for convolution. */
  void initialize(void)
  {
    Base::initialize();

    // TODO: Do we need this, alternatives?
    const double dt = TimeStep;

    // I'm assuming that these parameters just need to be of the
    // correct order of magnitude, but perhaps this is worth checking.
    const double eta = pow(mu0/eps0, 0.5);

    // To compute the step size at the first non stretched index
    //    TODO: This isn't safe
    unsigned dout = Forward ? Depth : Base::length(Axis) - Depth - 1;
    auto doutit = Base::next(Base::dataBegin(), Axis, dout);

    for (unsigned i = 0; i < Depth; ++i)
    {
      double dE, dH;
      if (Forward)
      {
        dE = (double) (Depth-i-0.5)/Depth;
        dH = (double) (Depth-i)/Depth;
      }
      else
      {
        dE = (double) (Depth-i)/Depth;
        dH = (double) (Depth-i-0.5)/Depth;
      }

      double smaxE = 1.0/(eta * Depth * Base::estep(doutit, Axis));
      double smaxH = 1.0/(eta * Depth * Base::hstep(doutit, Axis));

      double se = smaxE * Grading::sigma(dE);
      double sh = smaxH * Grading::sigma(dH);

      double ae = Grading::alpha(dE);
      double ah = Grading::alpha(dH);

      double ke = Grading::kappa(dE);
      double kh = Grading::kappa(dH);

      double Tol = 1.0e-6;

      m_eBeta[i] = exp(-(se/ke + ae)*(dt/eps0));
      m_hBeta[i] = exp(-(sh/kh + ah)*(dt/eps0));
      m_eAlpha[i] = se < Tol && ae < Tol ?
          m_eBeta[i]-1.0 : (m_eBeta[i]-1.0)*se/(se*ke + ke*ke*ae);
      m_hAlpha[i] = sh < Tol && ah < Tol ?
          m_hBeta[i]-1.0 : (m_hBeta[i]-1.0)*sh/(sh*kh + kh*kh*ah);
    }

    // Reset the convolution data
    for (auto& it : m_ephi) it = Utilities::Vec3d(0.0, 0.0, 0.0);
    for (auto& it : m_hphi) it = Utilities::Vec3d(0.0, 0.0, 0.0);
  }

  /** Update the E-field and the internal convolution data. */
  void advanceE(typename Base::TimingIterator it)
  {
    // To compute the step size at the first non stretched index
    //    TODO: This isn't safe
    unsigned dout = Forward ? Depth : Base::length(Axis) - Depth - 1;
    auto doutit = Base::next(Base::dataBegin(), Axis, dout);
    double hstep = Base::hstep(doutit, Axis);

    unsigned count = 0;

    for (unsigned i = imin; i < imax; ++i)
    {
      for (unsigned j = jmin; j < jmax; ++j)
      {
        for (unsigned k = kmin; k < kmax; ++k)
        {
          unsigned depth = Axis == 2 ? i : Axis == 1 ? j : k;
          if (!Forward) depth = Depth - (depth-amin) - 1;

          auto loc = Base::at(i, j, k);
          CpmlAxes::cpmlE(*this, it, depth, loc, count++, hstep);
        }
      }
    }

    Base::advanceE(it);
  }

  /** Update the H-field and the internal convolution data. */
  void advanceH(typename Base::TimingIterator it)
  {
    // To compute the step size at the first non stretched index
    //    TODO: This isn't safe
    unsigned dout = Forward ? Depth : Base::length(Axis) - Depth - 1;
    auto doutit = Base::next(Base::dataBegin(), Axis, dout);
    double estep = Base::estep(doutit, Axis);

    unsigned count = 0;

    for (unsigned i = imin; i < imax; ++i)
    {
      for (unsigned j = jmin; j < jmax; ++j)
      {
        for (unsigned k = kmin; k < kmax; ++k)
        {
          unsigned depth = Axis == 2 ? i : Axis == 1 ? j : k;
          if (!Forward) depth = Depth - (depth-amin) - 1;

          auto loc = Base::at(i, j, k);
          CpmlAxes::cpmlH(*this, it, depth, loc, count++, estep);
        }
      }
    }

    Base::advanceH(it);
  }
};

} // namespace Cpml
} // namespace Fdtd

#endif // #ifndef FDTD_CPML_CPML_HPP

