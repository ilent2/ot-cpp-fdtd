/* Output/ParameterSelectors.hpp - Common parameter selectors for output.
 *
 * Copyright (C) 2016 Isaac Lenton (aka ilent2)
 */

#ifndef FDTD_OUTPUT_PARAMETER_SELECTORS_HPP
#define FDTD_OUTPUT_PARAMETER_SELECTORS_HPP

#include <string>

namespace Fdtd {
namespace Output {
namespace Parameters {

/** A parameter selector for permittivity. */
struct Permittivity
{
  template <typename Base, typename It>
  static auto Parameter(Base& b, It loc)
      -> decltype(loc->permittivity())
  {
    return loc->permittivity();
  }

  static std::string Name(void) { return "Permittivity"; }
};

/** A parameter selector for permeability. */
struct Permeability
{
  template <typename Base, typename It>
  static auto Parameter(Base& b, It loc)
      -> decltype(loc->permebaility())
  {
    return loc->permeability();
  }

  static std::string Name(void) { return "Permeability"; }
};

/** A parameter selector for E-field data. */
struct Efield
{
  template <typename Base, typename It>
  static auto Parameter(Base& b, It loc)
      -> decltype(loc->efield())
  {
    return loc->efield();
  }

  static std::string Name(void) { return "E-field"; }
};

/** Select the E-field minus the time average E-field. */
struct EfieldMinusTimeAverage
{
  template <typename Base, typename It>
  static auto Parameter(Base& b, It loc)
      -> decltype(Base::acEfield(loc))
  {
    return Base::acEfield(loc);
  }

  static std::string Name(void) { return "E-field minus time average"; }
};

/** A parameter selector for H-field data. */
struct Hfield
{
  template <typename Base, typename It>
  static auto Parameter(Base& b, It loc)
      -> decltype(loc->hfield())
  {
    return loc->hfield();
  }

  static std::string Name(void) { return "H-field"; }
};

/** Select the Poynting vector. */
struct PoyntingMinusTimeAverage
{
  template <typename Base, typename It>
  static auto Parameter(Base& b, It loc)
      -> decltype(Base::acEfield(loc))
  {
    // TODO: This isn't right, the components are at different locations
    return Base::acEfield(loc).cross(Base::acHfield(loc));
  }

  static std::string Name(void) { return "Poynting minus time average"; }
};

} // namespace Parameters
} // namespace Output
} // namespace Fdtd

#endif // #ifndef FDTD_OUTPUT_PARAMTER_SELECTORS_HPP

