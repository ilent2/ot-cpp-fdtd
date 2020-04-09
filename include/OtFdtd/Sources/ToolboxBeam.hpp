/* Sources/ToolboxBeam.hpp - Optical Tweezers Toolbox beam.
 *
 * Loads files describing the point-matched coefficients computed
 * using the optical tweezers toolbox.
 *
 * The E and H fields are then calculated using the VSWF sum
 *
 *  E = \sum_n(1, N) \sum_m(-n, n) a_mn * M_mn + b_mn * N_mn
 *  E = -j * \sum_n(1, N) \sum_m(-n, n) b_mn * M_mn + a_mn * N_mn
 *
 * Copyright (C) 2016 Isaac Lenton (aka ilent2)
 */

#ifndef FDTD_SOURCES_TOOLBOX_BEAM_HPP
#define FDTD_SOURCES_TOOLBOX_BEAM_HPP

#include <string>
#include "Offset.hpp"
#include "Utilities/Vec.hpp"
#include <vector>
#include <fstream>

// From OTT toolbox
#include "emFieldXyz.hpp"

namespace Fdtd {
namespace Sources {

/** Load coefficients describing the contribution of each VSWF order.
 *
 * @tparam fnMn       : File name of file containing truncated m, n values.
 * @tparam fnAb       : File name of file complex amplitude coefficients.
 * @tparam Velocity   : Wave velocity in background medium.
 */
template <const std::string& fnMn, const std::string& fnAb,
    const double& Velocity, typename... CoordinateModifiers>
class ToolboxBeam
{
  std::vector<Utilities::Vec3d> m_elocs;
  std::vector<Utilities::Vec3d> m_hlocs;

  double m_dWavenumber;
  bool m_bCommonCoordinates;    ///< Use the same locations for E and H.

public:

  constexpr static bool Forward = true;

  /** Construct the beam object applying modifiers with Base to coordinates. */
  template <typename... Base>
  ToolboxBeam(const std::vector<Utilities::Vec3d>& elocs,
      const std::vector<Utilities::Vec3d>& hlocs, double AngularFrequency,
      const Base&... base)
  {
    static_assert(sizeof...(Base) <= 1, "Too many arguments");

    m_elocs.insert(m_elocs.end(), elocs.begin(), elocs.end());
    m_hlocs.insert(m_hlocs.end(), hlocs.begin(), hlocs.end());
    m_bCommonCoordinates = false;

    for (auto& it : m_elocs)
      it = Offset::Modify<CoordinateModifiers...>::modify(it, base...);
    for (auto& it : m_hlocs)
      it = Offset::Modify<CoordinateModifiers...>::modify(it, base...);

    m_dWavenumber = AngularFrequency / Velocity;
  }

  /** Construct the beam object applying modifiers with Base to coordinates. */
  template <typename... Base>
  ToolboxBeam(const std::vector<Utilities::Vec3d>& locs,
      double AngularFrequency, const Base&... base)
  {
    static_assert(sizeof...(Base) <= 1, "Too many arguments");

    m_elocs.insert(m_elocs.end(), locs.begin(), locs.end());
    m_bCommonCoordinates = true;

    for (auto& it : m_elocs)
      it = Offset::Modify<CoordinateModifiers...>::modify(it, base...);

    m_dWavenumber = AngularFrequency / Velocity;
  }

  /** Read the MN file. */
  std::vector<int> readMn(void) const
  {
    std::ifstream fp(fnMn);
    if (!fp)
    {
      throw std::runtime_error("Unable to open mn file: " + fnMn);
    }

    std::vector<int> ret;
    int v;
    while (fp >> v) ret.emplace_back(v);

    return ret;
  }

  /** Read the AB file. */
  static std::vector<std::complex<double>> readAb(void)
  {
    std::ifstream fp(fnAb);
    if (!fp)
    {
      throw std::runtime_error("Unable to open ab file: " + fnAb);
    }

    std::vector<std::complex<double>> ret;
    double re, im;
    while (fp >> re && fp >> im) ret.emplace_back(re, im);

    return ret;
  }

  /** Calculate the power of the beam from the beam coefficients. */
  static double calculatePower(double impedance=1.0, double wavenumber=1.0)
  {
    std::vector<std::complex<double>> ab = readAb();
    std::vector<std::complex<double>> a(ab.begin(), ab.begin()+ab.size()/2);
    std::vector<std::complex<double>> b(ab.begin()+ab.size()/2, ab.end());

    double ret = 0.0;
    for (unsigned i = 0; i < ab.size()/2; ++i)
    {
      ret += pow(std::abs(a[i]), 2.0) + pow(std::abs(b[i]), 2.0);
    }

    return ret/2.0/impedance/wavenumber/wavenumber;
  }

  void getFields(std::vector<Utilities::Vec3c>& efield,
      std::vector<Utilities::Vec3c>& hfield)
  {
    std::vector<int> mn = readMn();
    std::vector<int> n(mn.begin(), mn.begin()+mn.size()/2);
    std::vector<int> m(mn.begin()+mn.size()/2, mn.end());

    std::vector<std::complex<double>> ab = readAb();
    std::vector<std::complex<double>> a(ab.begin(), ab.begin()+ab.size()/2);
    std::vector<std::complex<double>> b(ab.begin()+ab.size()/2, ab.end());

    // Check the length is good
    if (mn.size() != ab.size())
    {
      throw std::runtime_error("Different lengths for ab and nm!");
    }

    // Offload to the Toolbox wrapper
    std::vector<Utilities::Vec3c> E;
    std::vector<Utilities::Vec3c> H;
    if (m_bCommonCoordinates)
    {
      std::tie(E, H) = Toolbox::emFieldXyz(m_elocs, m, n, a, b, m_dWavenumber);
    }
    else
    {
      E = Toolbox::emFieldXyzE(m_elocs, m, n, a, b, m_dWavenumber);
      H = Toolbox::emFieldXyzH(m_hlocs, m, n, a, b, m_dWavenumber);
    }

    // Unpack E/H
    efield.insert(efield.end(), E.begin(), E.end());
    hfield.insert(hfield.end(), H.begin(), H.end());
  }
};

} // namespace Sources
} /* namespace Fdtd */

#endif /* #ifndef FDTD_SOURCES_TOOLBOX_BEAM_HPP */

