/* Fdtd/Output/Report.hpp - Report information about the module stack.
 *
 * Purpose: To verbosely inform the user about the module configuration
 *    and when each simulation run starts and ends.
 *
 * Usage: The module implements three functions which are called
 *    on construction, after initialization and before finalization.
 *    To report information to the user during one of these phases
 *    simply implement the corresponding report_ method in any
 *    module class between (and excluding) the Fallback class
 *    and the EvaluateMethods class.
 *
 * Copyright (C) 2016 Isaac Lenton (aka ilent2)
 */

#ifndef FDTD_OUTPUT_REPORT_HPP
#define FDTD_OUTPUT_REPORT_HPP

#include "../EmptyModule.hpp"
#include <iostream>

namespace Fdtd {
namespace Output {

struct Report : public EmptyModule
{

  template <typename Base> struct Fallback : public Base
  {
    template <typename T> void report_construct(T& fp) const {}
    template <typename T> void report_initialize(T& fp) const {}
    template <typename T> void report_finalize(T& fp) const {}
  };

  template <typename Base> struct EvaluateMethods : public Base
  {
    /** Report construction. */
    EvaluateMethods(void)
    {
      std::cout << "Report: constructing..." << std::endl;
      Base::report_construct(std::cout);
      std::cout << "Report finished." << std::endl;
    }

    /** Report initialization. */
    void initialize(void)
    {
      Base::initialize();

      std::cout << "Report: starting simulation..." << std::endl;
      Base::report_initialize(std::cout);
      std::cout << "Report finished." << std::endl;
    }

    /** Report finalization. */
    void finalize(void)
    {
      std::cout << "Report: simulation finished..." << std::endl;
      Base::report_finalize(std::cout);
      std::cout << "Report finished." << std::endl;

      Base::finalize();
    }

  };

};

} // namespace Output
} // namespace Fdtd

#endif // #ifndef FDTD_OUTPUT_REPORT_HPP

