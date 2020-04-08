/* Fdtd/Materials/SimplePermeability.hpp - Simple isotropic permeability.
 *
 * Copyright (C) 2016 Isaac Lenton (aka ilent2)
 */

#ifndef FDTD_MATERIALS_SIMPLE_PERMEABILITY_HPP
#define FDTD_MATERIALS_SIMPLE_PERMEABILITY_HPP

#include "OneStepPermeability.hpp"

namespace Fdtd {
namespace Materials {

struct SimplePermeability : public _Permeability::OneStep
{

  template <typename Base> class ElementData
      : public _Permeability::OneStep::ElementData<Base>
  {
    double m_permeability;

  public:

    //@{
    /** Access the material permeability value. */
    const double& isoPermeability(void) const { return m_permeability; }
    const double& permeability(void) const { return m_permeability; }
    double& permeability(void) { return m_permeability; }
    //@}
  };

  template <typename Base> class EvolveMethods
      : public _Permeability::OneStep::EvolveMethods<Base>
  {
  public:

    /** Report our existence. */
    template <typename T> void report_construct(T& fp) const
    {
      fp << "SimplePermeability: " << std::endl;
      Base::report_construct(fp);
    }
  };

};

} // namespace Materials
} // namespace Fdtd

#endif // #ifndef FDTD_MATERIALS_SIMPLE_PERMEABILITY_HPP

