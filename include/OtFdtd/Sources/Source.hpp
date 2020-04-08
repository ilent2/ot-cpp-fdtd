/* Sources/Source.hpp - Implementations of pulsed and continuous wave source.
 *
 * These object take a beam object and calculate and store the
 * complex field values for the given locations.  Used by TFSF sources.
 *
 * Copyright (C) 2016 Isaac Lenton (aka ilent2)
 */

#include <vector>
#include <array>
#include "../Utilities/Vec.hpp"

#ifndef FDTD_SOURCES_SOURCE_HPP
#define FDTD_SOURCES_SOURCE_HPP

namespace Fdtd {
namespace Sources {

/** A continuous wave source.
 *
 * @tparam Beam       : Description of the beam.
 * @tparam Frequency  : Angular frequency of source.
 */
template <typename Beam, const double& Frequency>
class ContinuousWave
{
  std::vector<Utilities::Vec3c> m_efield;      ///< Complex E-field amplitudes
  std::vector<Utilities::Vec3c> m_hfield;      ///< Complex H-field amplitudes

public:

  /** Calculate the complex field amplitudes. */
  template <typename Base, typename LocType>
  void initialize(const Base& base,
      const LocType& elocs, const LocType& hlocs)
  {
    m_efield.clear();
    m_hfield.clear();
    Beam(elocs, hlocs, Frequency, base).getFields(m_efield, m_hfield);
  }

  /** Initialize using shared locations for E and H. */
  template <typename Base, typename LocType>
  void initialize(const Base& base, const LocType& locs)
  {
    m_efield.clear();
    m_hfield.clear();
    Beam(locs, Frequency, base).getFields(m_efield, m_hfield);
  }

  template <typename T> void report_initialize(T& fp) const
  {
    // TODO
  }

  /** Calculate and return the real field amplitudes at time. */
  std::vector<Utilities::Vec3d> efield(double time)
  {
    const double sign = Beam::Forward ? -1.0 : 1.0;
    std::complex<double> dt =
        exp(std::complex<double>(0.0, sign) * Frequency * time);

    std::vector<Utilities::Vec3d> ret(m_efield.size());
    for (unsigned i = 0; i < m_efield.size(); ++i)
    {
      ret[i] = real(dt * m_efield[i]);
    }

    return ret;
  }

  /** Calculate and return the real field amplitudes at time. */
  std::vector<Utilities::Vec3d> hfield(double time)
  {
    const double sign = Beam::Forward ? -1.0 : 1.0;
    std::complex<double> dt =
        exp(std::complex<double>(0.0, sign) * Frequency * time);

    std::vector<Utilities::Vec3d> ret(m_hfield.size());
    for (unsigned i = 0; i < m_hfield.size(); ++i)
    {
      ret[i] = real(dt * m_hfield[i]);
    }

    return ret;
  }

  /** Return a copy of the complex E-field data. */
  std::vector<Utilities::Vec3c> cefield(void) const
  {
    return m_efield;
  }

  /** Return a copy of the complex H-field data. */
  std::vector<Utilities::Vec3c> chfield(void) const
  {
    return m_hfield;
  }

};

/** A source containing multiple discrete frequency components.
 *
 * @tparam Beam         : Description of the incident beam.
 * @tparam Components   : Number of frequency components.
 * @tparam Frequencies  : Frequencies for each component.
 * @tparam Amplitudes   : Amplitude of each frequency component.
 */
template <typename Beam, unsigned Components,
    std::array<double, Components>& Frequencies,
    std::array<double, Components>& Amplitudes>
class MultiWave
{
  std::vector<Utilities::Vec3c> m_efield;  ///< E-field amplitudes
  std::vector<Utilities::Vec3c> m_hfield;  ///< H-field amplitudes
  unsigned m_unElocs;           ///< Number of E-Field locations
  unsigned m_unHlocs;           ///< Number of H-Field locations

public:
  /** Calculate the complex field amplitudes. */
  template <typename Base, typename LocType>
  void initialize(const Base& base,
      const LocType& elocs, const LocType& hlocs)
  {
    m_efield.clear();
    m_hfield.clear();
    m_unElocs = elocs.size();
    m_unHlocs = hlocs.size();

    for (unsigned i = 0; i < Components; ++i)
    {
      Beam(elocs, hlocs, Frequencies[i], base).getFields(m_efield, m_hfield);
    }
  }

  template <typename Base, typename LocType>
  void initialize(const Base& base, const LocType& locs)
  {
    m_efield.clear();
    m_hfield.clear();
    m_unElocs = locs.size();
    m_unHlocs = locs.size();

    for (unsigned i = 0; i < Components; ++i)
    {
      Beam(locs, Frequencies[i], base).getFields(m_efield, m_hfield);
    }
  }

  /** Calculate and return the real field amplitudes at time. */
  std::vector<Utilities::Vec3d> efield(double time)
  {
    auto efield = m_efield.cbegin();

    std::vector<Utilities::Vec3d> ret(m_unElocs, 0.0);
    for (unsigned i = 0; i < Components; ++i)
    {
      const double sign = Beam::Forward ? -1.0 : 1.0;
      std::complex<double> dt = Amplitudes[i] *
          exp(std::complex<double>(0.0, sign) * Frequencies[i] * time);

      for (unsigned j = 0; j < m_unElocs; ++j)
      {
        ret[j] += real(dt * *efield++);
      }
    }

    return ret;
  }

  /** Calculate and return the real field amplitudes at time. */
  std::vector<Utilities::Vec3d> hfield(double time)
  {
    auto hfield = m_hfield.cbegin();

    std::vector<Utilities::Vec3d> ret(m_unHlocs, 0.0);
    for (unsigned i = 0; i < Components; ++i)
    {
      const double sign = Beam::Forward ? -1.0 : 1.0;
      std::complex<double> dt = Amplitudes[i] *
          exp(std::complex<double>(0.0, sign) * Frequencies[i] * time);

      for (unsigned j = 0; j < m_unHlocs; ++j)
      {
        ret[j] += real(dt * *hfield++);
      }
    }

    return ret;
  }
};

/** Data for Pulsed. */
template <typename Beam, unsigned Components, const double& Spacing,
    const double& Central>
class _Pulsed
{
public:
  static std::array<double, Components> Frequencies;
  static std::array<double, Components> Amplitudes;

  /** Populate/update arrays with values calculated using Spacing & Central.
   *
   * We do this so Spacing and Central can be changed at runtime without
   * reallocation of the arrays.
   */
  static void initialize(void)
  {
    double pulseWidth = (Components-1.0) * Spacing;
    double lowerFrequency = Central - pulseWidth/2.0;
    for (unsigned i = 0; i < Components; ++i)
    {
      Frequencies[i] = lowerFrequency + i*Spacing;
      Amplitudes[i] = exp(-pow((Frequencies[i]-Central)/pulseWidth, 2.0));
    }
  }
};

/** GCC insists that these are defined. */
template <typename Beam, unsigned Components, const double& Spacing,
    const double& Central> std::array<double, Components>
    _Pulsed<Beam, Components, Spacing, Central>::Frequencies;

/** GCC insists that these are defined. */
template <typename Beam, unsigned Components, const double& Spacing,
    const double& Central> std::array<double, Components>
    _Pulsed<Beam, Components, Spacing, Central>::Amplitudes;

/** A Gaussian pulse source.
 *
 * This is a helper class for a MultiWave source containing the
 * discrete frequencies required to represent a Gaussian pulse.
 *
 *  @tparam Beam        : Object describing the source beam.
 *  @tparam Components  : Number of frequency components.
 *  @tparam Spacing     : Spacing between frequency components.
 *  @tparam Central     : Central frequency of pulse.
 */
template <typename Beam, unsigned Components, const double& Spacing,
    const double& Central>
class Pulsed : public MultiWave<Beam, Components,
    _Pulsed<Beam, Components, Spacing, Central>::Frequencies,
    _Pulsed<Beam, Components, Spacing, Central>::Amplitudes>
{
  using BaseType = MultiWave<Beam, Components,
      _Pulsed<Beam, Components, Spacing, Central>::Frequencies,
      _Pulsed<Beam, Components, Spacing, Central>::Amplitudes>;

public:
  template <typename Base, typename LocType>
  void initialize(const Base& base,
      const LocType& elocs, const LocType& hlocs)
  {
    _Pulsed<Beam, Components, Spacing, Central>::initialize();
    BaseType::initialize(base, elocs, hlocs);
  }

  template <typename T> void report_initialize(T& fp) const
  {
    fp << "Using Pulsed beam with " << Components
        << " components." << std::endl;
  }
};

} // namespace Sources
} // namespace Fdtd

#endif // #ifndef FDTD_SOURCES_SOURCE_HPP

