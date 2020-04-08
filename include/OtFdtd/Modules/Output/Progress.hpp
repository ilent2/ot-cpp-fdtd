/* Output/Progress.hpp - Output the progress every N cycles.
 *
 * Copyright (C) 2016 Isaac Lenton (aka ilent2)
 */

#ifndef FDTD_OUTPUT_PROGRESS_HPP
#define FDTD_OUTPUT_PROGRESS_HPP

#include "Modules/EmptyModule.hpp"

namespace Fdtd {
namespace Output {

/** Output the progress every N output cycles.
 *
 *  @tparam StrideWindow    : The output stride description.
 *  @tparam FileName        : Filename to write to.
 */
template <unsigned N>
struct Progress : public Fdtd::EmptyModule
{
  template <typename Base> class EvolveMethods : public Base
  {
  public:

    /** Report our existence in construction. */
    template <typename T> void report_initialize(T& fp) const
    {
      fp << "Output:Progress outputting every " << N << std::endl;
      Base::report_initialize(fp);
    }

    /** Write the time to file. */
    void inspectE(typename Base::TimingIterator it)
    {
      if (it->index() % N == 0)
      {
        std::cout << "Progress... " << it->index() << std::endl;
      }

      // Continue inspection
      Base::inspectE(it);
    }
  };
};

} // namespace Output
} // namespace Fdtd

#endif // #ifndef FDTD_OUTPUT_PROGRESS_HPP

