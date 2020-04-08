/* Fdtd/Materials/SimplePermittivity.hpp - Simple isotropic permittivity.
 *
 * Copyright (C) 2016 Isaac Lenton (aka ilent2)
 */

#ifndef FDTD_MATERIALS_SIMPLE_PERMITTIVITY_HPP
#define FDTD_MATERIALS_SIMPLE_PERMITTIVITY_HPP

#include "OneStepPermittivity.hpp"

namespace Fdtd {
namespace Materials {

struct SimplePermittivity : public _Permittivity::OneStep
{

  template <typename Base> class ElementData
      : public _Permittivity::OneStep::ElementData<Base>
  {
    double m_permittivity;

  public:

    //@{
    /** Access the material permittivity value. */
    const double& isoPermittivity(void) const { return m_permittivity; }
    const double& permittivity(void) const { return m_permittivity; }
    double& permittivity(void) { return m_permittivity; }
    //@}
  };

  template <typename Base> class EvolveMethods
      : public _Permittivity::OneStep::EvolveMethods<Base>
  {
  public:

    /** Report our existence. */
    template <typename T> void report_construct(T& fp) const
    {
      fp << "Using SimplePermittivity." << std::endl;
      Base::report_construct(fp);
    }
  };

};


/*
  // TODO: Do we need to use AdvanceMetods::addDfield still?

  template <typename Base>
  class AdvanceMethods : public Base
  {
  public:

    template <typename T>
    void addDfield(unsigned e, T value)
    {
      auto& elm = Base::atIndex(e);
      elm.efield() += value / elm.permittivity();
    }

  };

*/

} // namespace Materials
} // namespace Fdtd

#endif // #ifndef FDTD_MATERIALS_SIMPLE_PERMITTIVITY_HPP

