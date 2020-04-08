/** Sources/PlaneWaveBeam.hpp - A simple plane wave beam.
 *
 * Copyright (C) 2016 Isaac Lenton (aka ilent2)
 */

#ifndef FDTD_SOURCES_PLANE_WAVE_BEAM_HPP
#define FDTD_SOURCES_PLANE_WAVE_BEAM_HPP

#include <vector>
#include "Offset.hpp"
#include "../Utilities/Vec.hpp"

namespace Fdtd {
namespace Sources {

/** Linearly polarized plane wave beam.
 *
 * @tparam Normal       : Direction of propagation.
 * @tparam Polarization : Direction of E-field.
 * @tparam Intensity    : Beam intensity.
 * @tparam Velocity     : Wave speed in medium.
 * @tparam CoordinateModifiers... : Modifiers to apply to the coordinates.
 */
template <const Utilities::Vec3d& Normal, const Utilities::Vec3d& Polarization,
    const double& Intensity, const double& Velocity,
    typename... CoordinateModifiers>
class PlaneWaveBeam
{
protected:
  std::vector<Utilities::Vec3d> locs;
  unsigned sizeElocs;

  Utilities::Vec3d m_normal;     ///< Propagation direction
  Utilities::Vec3c m_enorm;      ///< E-Field direction
  Utilities::Vec3c m_hnorm;      ///< H-Field direction

  double m_dWavenumber;   ///< Wavenumber

public:

  /** Always true, if we want to flip the direction, flip the B field. */
  constexpr static bool Forward = true;

  /** Construct the beam object applying modifiers to coordinates. */
  template <typename... Base>
  PlaneWaveBeam(const std::vector<Utilities::Vec3d>& elocs,
      const std::vector<Utilities::Vec3d>& hlocs, double AngularFrequency,
      const Base&... base)
  {
    sizeElocs = elocs.size();
    locs.insert(locs.end(), elocs.begin(), elocs.end());
    locs.insert(locs.end(), hlocs.begin(), hlocs.end());

    m_dWavenumber = AngularFrequency / Velocity;

    for (auto& it : locs)
      it = Offset::Modify<CoordinateModifiers...>::modify(it, base...);

    // Compute the normalized normal
    m_normal = Normal.normalized();
    m_hnorm = m_normal.cross(Polarization.normalized());
    m_enorm = -m_hnorm.cross((Utilities::Vec3c) m_normal);
  }

  std::complex<double> compute(unsigned i)
  {
    double dz = locs[i].dot(m_normal);
    return Intensity * exp(std::complex<double>(0.0, 1.0)*m_dWavenumber*dz);
  }

  void getFields(std::vector<Utilities::Vec3c>& efield,
      std::vector<Utilities::Vec3c>& hfield)
  {
    for (unsigned i = 0; i < sizeElocs; ++i)
    {
      efield.push_back(compute(i) * m_enorm);
    }
    for (unsigned i = sizeElocs; i < locs.size(); ++i)
    {
      hfield.push_back(compute(i) * m_hnorm);
    }
  }
};

/** Elliptically polarized plane wave beam.
 *
 * @tparam Normal       : Direction of propagation.
 * @tparam Polarization : Primary direction of E-field.
 * @tparam IntensityP   : Primary beam intensity.
 * @tparam IntensityS   : Secondary beam intensity.
 * @tparam Velocity     : Wave speed in medium.
 * @tparam CoordinateModifiers... : Modifiers to apply to the coordinates.
 */
template <const Utilities::Vec3d& Normal, const Utilities::Vec3d& Polarization,
    const double& IntensityP, const double& IntensityS, const double& Velocity,
    typename... CoordinateModifiers>
class EllipticalPlaneWaveBeam
  : public PlaneWaveBeam<Normal, Polarization,
      IntensityP, Velocity, CoordinateModifiers...>
{
  using Beam = PlaneWaveBeam<Normal, Polarization,
      IntensityP, Velocity, CoordinateModifiers...>;

public:
  /** Construct the beam object applying modifiers to coordinates. */
  template <typename... Base>
  EllipticalPlaneWaveBeam(const std::vector<Utilities::Vec3d>& elocs,
      const std::vector<Utilities::Vec3d>& hlocs, double AngularFrequency,
      const Base&... base)
    : Beam(elocs, hlocs, AngularFrequency, base...)
  {
  }

  void getFields(std::vector<Utilities::Vec3c>& efield,
      std::vector<Utilities::Vec3c>& hfield)
  {
    std::complex<double> shift
        = exp(std::complex<double>(0.0, M_PI/2.0)) * IntensityS/IntensityP;

    for (unsigned i = 0; i < Beam::sizeElocs; ++i)
    {
      efield.push_back(Beam::compute(i) * Beam::m_enorm
          + shift * Beam::compute(i) * Beam::m_hnorm);
    }
    for (unsigned i = Beam::sizeElocs; i < Beam::locs.size(); ++i)
    {
      hfield.push_back(Beam::compute(i) * Beam::m_hnorm
          - shift * Beam::compute(i) * Beam::m_enorm);
    }
  }
};

/** Circularly polarized plane wave beam.
 * @tparam Normal       : Direction of propagation.
 * @tparam Polarization : Initial direction of E-field.
 * @tparam Intensity    : Amplitude for the E-field.
 * @tparam Velocity     : Wave speed in medium.
 * @tparam CoordinateModifiers... : Modifiers to apply to the coordinates.
 */
template <const Utilities::Vec3d& Normal, const Utilities::Vec3d& Polarization,
    const double& Intensity, const double& Velocity,
    typename... CoordinateModifiers>
using CircularPlaneWaveBeam = EllipticalPlaneWaveBeam<Normal,
    Polarization, Intensity, Intensity, Velocity, CoordinateModifiers...>;

} // namespace Sources
} // namespace Fdtd

#endif /* #ifndef FDTD_SOURCES_PLANE_WAVE_BEAM_HPP */

