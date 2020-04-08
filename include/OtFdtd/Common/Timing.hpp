/* Common/Timing.hpp - Types for commonly used simulation time resolutions.
 *
 * Copyright (C) 2016 Isaac Lenton (aka ilent2)
 */

#ifndef FDTD_COMMON_TIMING_HPP
#define FDTD_COMMON_TIMING_HPP

#include "Modules/Timing/Linear.hpp"

namespace Common {

/** Compute the time step size and integer cycle length.
 *
 * @tparam material       : An object with a frequency member constant.
 * @tparam max_time_step  : The maximum time step size permissible.
 * @tparam multiple       : Cycle length should be divisible by this number.
 */
template <const double& max_time_step, typename material, unsigned Multiple=3>
class DivisibleTimeStep
{
  static constexpr double min_cycle_length
      = 1.0 / max_time_step / material::frequency;

public:

  /** Expose the multiple as an interesting parameter. */
  static constexpr unsigned multiple = Multiple;

  /** The computed integer optical cycle length. */
  static constexpr unsigned cycle_length
      = multiple * ceil(min_cycle_length / multiple);

  /** The computed time step size. */
  static constexpr double step_size
      = 1.0 / cycle_length / material::frequency;

  /** A useful parameter for ForceTorque strides. */
  static constexpr unsigned output_stride = cycle_length / multiple;
};

namespace Timing {

/** A helper providing Linear time steps with a multiple. */
template <unsigned num_cycles, const double& max_step_size,
    typename material, unsigned multiple=3>
class Linear : public DivisibleTimeStep<max_step_size, material, multiple>,
    public Fdtd::Timing::Linear<num_cycles*
        DivisibleTimeStep<max_step_size, material, multiple>::cycle_length,
        DivisibleTimeStep<max_step_size, material, multiple>::step_size>
{
public:
  static constexpr unsigned total_steps = num_cycles *
      DivisibleTimeStep<max_step_size, material, multiple>::cycle_length;
};

} // namespace Timing
} // namespace Common

#endif // #ifndef FDTD_COMMON_TIMING_HPP

