/* StrideWindow.hpp - Description of output stride, window, start and end.
 *
 * Copyright (C) 2016 Isaac Lenton (aka ilent2)
 */

#ifndef FDTD_OUTPUT_STRIDE_WINDOW_HPP
#define FDTD_OUTPUT_STRIDE_WINDOW_HPP

namespace Fdtd {
namespace Output {

/** Description of the output stride, window, start and end.
 *
 * @tparam Stride     : Time steps between each calculation.
 * @tparam Window     : Number of time steps to average over.
 * @tparam Start      : Time step to start calculation.
 * @tparam End        : Time step to end calculation.
 *
 * For force/torque calculations, Window*Stride will typically be a
 * multiple of the optical frequency.
 */
template <unsigned stride=1, unsigned window=1, int start=-1, int end=-1>
struct StrideWindow
{
  static constexpr unsigned Stride = stride;
  static constexpr unsigned Window = window;
  static constexpr int Start = start;
  static constexpr int End = end;

  /** Determine if index is an output index.
   *
   * Returns true for the first Window indices before Start.
   */
  static bool outputNow(unsigned index)
  {
    // TODO: This might actually be a bit slow...
    //    Up to 50% total runtime in one test (might be problem elsewhere)
    //    Perhaps we should look at alternatives?

    if (Start >= 0 && (int) index <= Start-int (Window*Stride)) return false;
    if (End > 0 && (int) index >= End) return false;
    if (index % Stride != 0) return false;
    return true;
  }

  /** Determine if index is an output index (with no Window).
   *
   * Returns true for the first Window indices before Start.
   */
  static bool outputNowNoWindow(unsigned index)
  {
    // TODO: This might actually be a bit slow...
    //    Up to 50% total runtime in one test (might be problem elsewhere)
    //    Perhaps we should look at alternatives?

    if (Start >= 0 && (int) index < Start) return false;
    if (End > 0 && (int) index >= End) return false;
    if (index % Stride != 0) return false;
    return true;
  }

  /** Helper for types with an index method (such as timing iterator). */
  template <typename It>
  static bool outputNow(const It& it)
  {
    return outputNow(it->index());
  }

  /** Helper for types with an index method (such as timing iterator). */
  template <typename It>
  static bool outputNowNoWindow(const It& it)
  {
    return outputNowNoWindow(it->index());
  }

  /** Create a string representation of this object. */
  static std::string toString(void)
  {
    std::stringstream ss;
    ss << "(" << Stride << ", " << Window
        << ", " << Start << ", " << End << ")";
    return ss.str();
  }
};

/** A helper for stride only stride windows. */
template <unsigned stride=1, int start=-1, int end=-1>
using StrideOnly = StrideWindow<stride, 1, start, end>;

} // namespace Output
} // namespace Fdtd

#endif // #ifndef FDTD_OUTPUT_STRIDE_HPP

