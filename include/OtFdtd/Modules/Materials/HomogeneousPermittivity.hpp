/* Fdtd/Materials/HomogeneousPermittivity.hpp - Homogeneous permittivity.
 *
 * Purpose: Provide a simple 1-step E-field update method with
 *    a homogeneous scalar permittivity.
 *
 * Copyright (C) 2016 Isaac Lenton (aka ilent2)
 */

#ifndef FDTD_MATERIALS_HOMOGENEOUS_PERMITTIVITY_HPP
#define FDTD_MATERIALS_HOMOGENEOUS_PERMITTIVITY_HPP

#include "OneStepPermittivity.hpp"

namespace Fdtd {
namespace Materials {

template <const double& Value>
struct HomogeneousPermittivity : public _Permittivity::OneStep
{

  template <typename Base> class ElementData
      : public _Permittivity::OneStep::ElementData<Base>
  {
  public:

    //@{
    /** Provide a method for accessing the permittivity. */
    const double& permittivity(void) const { return Value; }
    const double& isoPermittivity(void) const { return Value; }
    //@}
  };

  template <typename Base> class EvolveMethods
      : public _Permittivity::OneStep::EvolveMethods<Base>
  {
  public:

    /** Report our existence. */
    template <typename T> void report_initialize(T& fp) const
    {
      fp << "HomogeneousPermittivity: value = " << Value << std::endl;
      Base::report_initialize(fp);
    }
  };

};

} // namespace Materials
} // namespace Fdtd

#endif // #ifndef FDTD_HOMOGENEOUS_PERMITTIVITY_HPP

