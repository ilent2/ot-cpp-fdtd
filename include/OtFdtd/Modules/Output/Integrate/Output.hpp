/* Output/Integrate/Output.hpp - Output targets for Integrate method.
 *
 * Copyright (C) 2016 Isaac Lenton (aka ilent2)
 */

#ifndef FDTD_OUTPUT_INTEGRATE_OUTPUT_HPP
#define FDTD_OUTPUT_INTEGRATE_OUTPUT_HPP

#include <fstream>

namespace Fdtd {
namespace Output {
namespace Integrate {
namespace Output {

/** Output for every value of the simulation. */
template <const std::string& FileName>
class FileEverything
{
  std::ofstream m_output;       // The output file

public:

  FileEverything(void)
  {
    // Nothing to do
  }

  void initialize(void)
  {
    // Try to open the output file
    m_output.open(FileName);
    if (!m_output)
      throw std::logic_error("Unable to open output file: " + FileName);

    // Write the file header
    m_output.precision(16);
    m_output << "# Index    Time     Value" << std::endl;
  }

  void finalize(void)
  {
    // Close the output file
    m_output.close();
  }

  /** Write the value to file. */
  template <typename Time, typename Value>
  void report(Time it, Value value)
  {
    m_output << it->index() << '\t'
        << it->efieldT() << '\t'
        << value << std::endl;
  }
};

/** Output the last value of every simulation.
 *
 * @tparam FileName     : The output file name.
 * @tparam T            : The type for the index variable.
 * @tparam Value        : The value for the index variable.
 * @tparam Increment    : Increment the index variable.
 * @tparam Append       : Append to the file instead of replacing.
 * @tparam AutoOpen     : Automatically open the file on construction.
 */
template <const std::string& FileName, typename T,
    const T& Value, bool Increment=false, bool Append=false,
    bool AutoOpen=true>
class FileLast
{
  static std::ofstream m_output;       // The output file

  unsigned m_ilast;
  double m_tlast;
  Utilities::Vec3d m_last;

public:

  /** Open the file. */
  static void open_file(bool append=Append)
  {
    // Try to open the output file
    m_output.close();
    m_output.open(FileName, append ? std::ios::app : std::ios::out);
    if (!m_output)
      throw std::logic_error("Unable to open output file: " + FileName);

    // Write the file header
    m_output.precision(16);
    if (!append) m_output << "# Parameter\tIndex\tTime\tValue" << std::endl;
  }

  FileLast(void)
  {
    if (AutoOpen) open_file();
  }

  void initialize(void)
  {
    // Nothing to do
  }

  void finalize(void)
  {
    // Write the last value
    m_output << Value << '\t'
        << m_ilast << '\t' << m_tlast << '\t' << m_last << std::endl;
  }

  /** Write the value to file. */
  template <typename Time, typename ValueT>
  void report(Time it, ValueT value)
  {
    m_last = value;
    m_ilast = it->index();
    m_tlast = it->efieldT();
  }
};

// Declare the output file for FileLast
template <const std::string& FileName, typename T,
    const T& Value, bool Increment, bool Append, bool AutoOpen>
std::ofstream FileLast<FileName, T, Value, Increment,
    Append, AutoOpen>::m_output;

/** Store the last value reported. */
template <typename T, T& Value>
struct StoreLast
{
  void initialize(void) {}
  void finalize(void) {}

  template <typename Time>
  void report(Time it, T value)
  {
    Value = value;
  }
};

/** An adapter to normalize the value. */
template <typename Base, typename T, const T& Value>
struct Normalized : public Base
{
  template <typename Time, typename ValueT>
  void report(Time it, ValueT value)
  {
    Base::report(it, value / Value);
  }
};

} // namespace Output (inside Integrate)
} // namespace Integrate
} // namespace Output
} // namespace Fdtd

#endif // #ifndef FDTD_OUTPUT_INTEGRATE_OUTPUT_HPP

