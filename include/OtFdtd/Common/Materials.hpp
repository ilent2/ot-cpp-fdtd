/* Common/Materials.hpp - Common material definitions.
 *
 * Copyright (C) 2016 Isaac Lenton (aka ilent2)
 */

#ifndef FDTD_COMMON_MATERIALS_HPP
#define FDTD_COMMON_MATERIALS_HPP

#include "Utilities/MinMaxAttribute.hpp"

namespace Common {

/** Common wavelengths. */
struct Wavelengths
{
  static constexpr double si1064 = 1064.0e-9;         //< (m)
};

/** A Vacuum material with SI units. */
template <const double& Wavelength = Wavelengths::si1064>
struct SiVacuum
{
  static constexpr double permittivity = 8.85e-12;    //< (F/m)
  static constexpr double permeability = 1.26e-6;     //< (H/m)
  static constexpr double impedance = pow(permittivity/permeability, 0.5);
  static constexpr double wavelength = Wavelength;
  static constexpr double wavenumber = 2.0 * M_PI / wavelength;
  static constexpr double speed = pow(permittivity*permeability, -0.5);
  static constexpr double frequency = speed / wavelength;
  static constexpr double angular_frequency = 2.0 * M_PI * frequency;
};

/** A Vacuum material with unitary permittivity/permeability/wavelength. */
struct UnitVacuum
{
  static constexpr double permittivity = 1.0;    //< (F/m)
  static constexpr double permeability = 1.0;     //< (H/m)
  static constexpr double impedance = pow(permittivity/permeability, -0.5);
  static constexpr double wavelength = 1.0;
  static constexpr double wavenumber = 2.0 * M_PI / wavelength;
  static constexpr double speed = pow(permittivity*permeability, -0.5);
  static constexpr double frequency = speed / wavelength;
  static constexpr double angular_frequency = 2.0 * M_PI * frequency;
};

/** A dielectric material specified by a refractive index. */
template <const double& Index, typename Vacuum = SiVacuum<>>
struct Dielectric
{
  static constexpr double index = Index;
  static constexpr double permittivity = Vacuum::permittivity*pow(index, 2.0);
  static constexpr double permeability = Vacuum::permeability;
  static constexpr double impedance = pow(permittivity/permeability, -0.5);
  static constexpr double wavelength = Vacuum::wavelength/index;
  static constexpr double wavenumber = 2.0 * M_PI / wavelength;
  static constexpr double speed = Vacuum::speed/index;
  static constexpr double frequency = Vacuum::frequency;
  static constexpr double angular_frequency = Vacuum::angular_frequency;
};

/** Common refractive indexes. */
struct Indices
{
  static constexpr double water = 1.33;
  static constexpr double polystyrene = 1.59;
  static constexpr double glass = 2.0;

  // Refractive indices for vaterite unit cell (from Wikipedia)
  static constexpr double vaterite_o = 1.55;
  static constexpr double vaterite_e = 1.65;
};

/** Attribute for wavelength. */
template <typename T>
struct AttrWavelength { static constexpr double value = T::wavelength; };

/** Calculate the minimum wavelength. */
template <typename... Types>
using wavelength_min_max = Utilities::MinMaxAttribute<
    AttrWavelength, Types...>;

/** Attribute for speed. */
template <typename T>
struct AttrSpeed { static constexpr double value = T::speed; };

/** Calculate the minimum wavelength. */
template <typename... Types>
using speed_min_max = Utilities::MinMaxAttribute<AttrSpeed, Types...>;

/** Estimate the maximum time step using the CFL condition.
 *
 * Assumes the spacing is used for all materials.
 */
template <const double& spacing, unsigned dimensions, typename... Types>
struct MaxTimeStepSize
{
  static constexpr double value(void)
  {
    return spacing/(speed_min_max<Types...>::max
        * pow((double) dimensions, 0.5));
  }
};

} // namespace Common

#endif // #ifndef FDTD_COMMON_MATERIALS_HPP

