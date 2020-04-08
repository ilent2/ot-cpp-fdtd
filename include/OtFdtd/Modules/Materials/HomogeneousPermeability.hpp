/* Fdtd/Materials/HomogeneousPermeability.hpp - Homogeneous permeability.
 *
 * Purpose: Provide a simple 1-step H-field update method with
 *    a homogeneous scalar permeability.
 *
 * Copyright (C) 2016 Isaac Lenton (aka ilent2)
 */

#ifndef FDTD_MATERIALS_HOMOGENEOUS_PERMEABILITY_HPP
#define FDTD_MATERIALS_HOMOGENEOUS_PERMEABILITY_HPP

#include "OneStepPermeability.hpp"

namespace Fdtd {
namespace Materials {

template <const double& Value>
struct HomogeneousPermeability : public _Permeability::OneStep
{

  template <typename Base> class ElementData
      : public _Permeability::OneStep::ElementData<Base>
  {
  public:

    //@{
    /** Provide a method for accessing the permeability. */
    const double& permeability(void) const { return Value; }
    const double& isoPermeability(void) const { return Value; }
    //@}
  };

  template <typename Base> class EvolveMethods
      : public _Permeability::OneStep::EvolveMethods<Base>
  {
  public:

    /** Report our existence. */
    template <typename T> void report_initialize(T& fp) const
    {
      fp << "HomogeneousPermeability: value = " << Value << std::endl;
      Base::report_initialize(fp);
    }
  };

};

} // namespace Materials
} // namespace Fdtd

#endif // #ifndef FDTD_HOMOGENEOUS_PERMEABILITY_HPP

