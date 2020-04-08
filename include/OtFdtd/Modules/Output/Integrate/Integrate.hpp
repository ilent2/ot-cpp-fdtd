/* Output/Integrate/Integrate.hpp - Integration methods, ready for use.
 *
 * See also:
 *    Integrate/Output.hpp
 *    Integrate/Parameters.hpp
 *    Integrate/Iterators.hpp
 *    StrideWindow.hpp
 *
 * Copyright (C) 2016 Isaac Lenton (aka ilent2)
 */

#ifndef FDTD_OUTPUT_INTEGRATE_INTEGRATE_HPP
#define FDTD_OUTPUT_INTEGRATE_INTEGRATE_HPP

#include <list>
#include "Modules/EmptyModule.hpp"

namespace Fdtd {
namespace Output {
namespace Integrate {

/** Integrate a parameter over a surface/region.
 *
 * @tparam Output       : An object describing the output destination.
 * @tparam Parameter    : The parameter to integrate.
 * @tparam Iterator     : A type describing the integration surface/volume.
 * @tparam Stride       : A type describing the output/calculation times.
 * @tparam Inject       : Other modules to inject before us.
 *
 * TODO: Use Inject to specify SubtractDc and other nice stuff
 */
template <typename Output, typename Parameter,
    typename Iterator, typename Stride, typename... Inject>
struct Integrate : public EmptyModule
{
  template <typename B> class EvolveMethods
      : public Iterator::template Iterator<B, Parameter>
  {
    using Base = typename Iterator::template Iterator<B, Parameter>;

		Utilities::Vec3d m_value;             ///< Moving window sum
    std::list<Utilities::Vec3d> m_window; ///< Moving window values
    Output m_output;                      ///< The output destination.

    void outputNow(typename Base::TimingIterator it)
    {
      // Do the integration at this time step
      auto value = Base::Integrate(it);

      // Add the values to the moving window
      m_window.push_back(value);
      m_value += m_window.back();

      // Ensure the moving window is not too long
      if (m_window.size() > Stride::Window)
      {
        m_value -= m_window.front();
        m_window.pop_front();
      }

      // Output the value from the moving window
      if ((int) it->index() >= Stride::Start)
      {
        m_output.report(it, m_value / Stride::Window);
      }
    }

  public:

    /** Initialize the moving window and underlying types. */
		void initialize(void)
    {
      Base::initialize();

      m_value = Utilities::Vec3d(0.0, 0.0, 0.0);
      m_window.clear();

      m_output.initialize();
    }

    /** Finalize the underlying types. */
    void finalize(void)
    {
      m_output.finalize();
      Base::finalize();
    }

    /** Report our existence in construction. */
    template <typename T> void report_construct(T& fp) const
    {
      fp << "Output:Integrate:"
          << " with stride " << Stride::toString()
          << " with iterator " << Iterator::repr()
          << " with kernel " << Parameter::repr() << std::endl;
      Base::report_construct(fp);
    }

    /** Integrate after the E field has been updated. */
    void inspectE(typename Base::TimingIterator it)
    {
      // Output if we are in an output phase
      if (Stride::outputNow(it)) outputNow(it);

      // Continue inspection
      Base::inspectE(it);
    }
  };
};

} // namespace Integrate
} // namespace Output
} // namespace Fdtd

#endif // #ifndef FDTD_OUTPUT_INTEGRATE_INTEGRATE_HPP

