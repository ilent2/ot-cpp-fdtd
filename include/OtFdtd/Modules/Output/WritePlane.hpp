/* Output/WritePlane.hpp - Write a 2-D slice of a field.
 *
 * Copyright (C) 2016 Isaac Lenton (aka ilent2)
 */

#ifndef FDTD_OUTPUT_WRITE_PLANE_HPP
#define FDTD_OUTPUT_WRITE_PLANE_HPP

#include "../EmptyModule.hpp"
#include <string>
#include <fstream>

namespace Fdtd {
namespace Output {

/** Output a 2-D slice of the specified field parameter to file.
 *
 * @tparam OutputFile : A output file with a %d substitution token.
 * @tparam Parameter  : The field parameter to output.
 * @tparam Axis       : Axis for output.
 * @tparam Depth      : Depth of plane along axis.
 * @tparam Stride     : Description of output range.
 *
 * The parameter type should have a single static member
 * taking two arguments: the assembled base class instance,
 * and the current location iterator.
 * The function should return the value at the specified element.
 * The parameter type should also have a streamable Name member.
 */
template <const std::string& OutputFile, typename Parameter,
    int Axis, unsigned Depth, typename Stride>
class WritePlane : public EmptyModule
{
public:
  template <typename Base> class EvolveMethods : public Base
  {
    static_assert(Depth < Base::length(Axis),
        "WritePlane: Depth must be less than axis length.");

    /** A recursive class to write the data. */
    template <int ThisAxis, typename... Args>
    class Writer
    {
    public:
      Writer(Base& base, std::ofstream& fp, Args... args)
      {
        WriteNow(base, fp, args...);
      }

      void WriteIndices(std::ofstream& fp)
      {
        // Terminal case: nothing to do
      }

      template <typename T, typename... Other>
      void WriteIndices(std::ofstream& fp, T first, Other... args)
      {
        fp << first << '\t';
        WriteIndices(fp, args...);
      }

      template <int A = ThisAxis>
      typename std::enable_if_t<(A > 0) && A-1 != Axis>
      WriteNow(Base& base, std::ofstream& fp, Args... args)
      {
        for (unsigned i = 0; i < base.length(A-1); ++i)
        {
          Writer<A-1, Args..., unsigned>(base, fp, args..., i);
        }
        fp << '\n';
      }

      template <int A = ThisAxis>
      typename std::enable_if_t<(A > 0) && A-1 == Axis>
      WriteNow(Base& base, std::ofstream& fp, Args... args)
      {
        Writer<A-1, Args..., unsigned>(base, fp, args..., Depth);
      }

      template <int A = ThisAxis>
      typename std::enable_if_t<A == 0>
      WriteNow(Base& base, std::ofstream& fp, Args... args)
      {
        WriteIndices(fp, args...);
        // TODO: Would it also be useful to write the coordinate?
        fp << Parameter::Parameter(base, base.at(args...)) << '\n';
      }
    };

  public:

    /** Report our existence. */
    template <typename T>
    void report_construct(T& fp) const
    {
      fp << "Output:WritePlane: Outputting " << Parameter::Name()
          << " to " << OutputFile
          << " on axis " << Axis << " at depth " << Depth
          << " with stride " << Stride::toString() << std::endl;
      Base::report_construct(fp);
    }

    /** Output one the E-field has been updated (end of iteration). */
    void inspectE(typename Base::TimingIterator it)
    {
      // We don't modify any state, so let the children do their stuff first
      Base::inspectE(it);

      if (Stride::outputNowNoWindow(it))
      {
        // Do string substitution and open file
        char buf[128];
        snprintf(buf, 128, OutputFile.c_str(), it->index());
        std::ofstream fp(buf, std::ios::binary);

        // Check we opened the file ok
        if (!fp)
        {
          throw std::runtime_error("Unable to open file: " + std::string(buf));
        }

        // Begin recursively iterating over axes/elements
        Writer<Base::AxisIndex>(*this, fp);
      }
    }
  };
};

} // namespace Output
} // namespace Fdtd

#endif // #ifndef FDTD_OUTPUT_WRITE_PLANE_HPP

