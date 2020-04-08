/* Output/Integrate/MaxwellStressTensor.hpp - Maxwell stress tensor.
 *
 * Copyright (C) 2016 Isaac Lenton (aka ilent2)
 */

#ifndef FDTD_OUTPUT_INTEGRATE_MAXWELL_STRESS_TENSOR_HPP
#define FDTD_OUTPUT_INTEGRATE_MAXWELL_STRESS_TENSOR_HPP

#include "StressTensorCommon.hpp"

namespace Fdtd {
namespace Output {
namespace Integrate {

/** Implementation of the Maxwell stress tensor for calculating Force/Torque.
 *
 * @tparam SubtractDc     : Subtract the DC field before calculation.
 */
template <bool SubtractDc=true>
class MaxwellStressTensor : public StressTensorCommon<SubtractDc>
{

  using MBase = StressTensorCommon<SubtractDc>;

  /** Calculate the Maxwell Stress tensor for vectors avgE and avgH.
   *
   * Compatible with anisotropic permittivity/permeability.
   */
  struct StressTensor
  {
    template <typename Eps, typename Mu>
    static Utilities::Mat3d stressTensor(
        Utilities::Vec3d avgE, Utilities::Vec3d avgH,
        Eps eps, Mu mu)
    {
      auto avgD = eps * avgE;
      auto avgB = mu * avgH;

      auto E2 = avgD.dot(avgE);
      auto H2 = avgB.dot(avgH);

      Utilities::Mat3d ret;
      ret.xx = (avgD.x * avgE.x - 0.5 * E2)
             + (avgB.x * avgH.x - 0.5 * H2);
      ret.yy = (avgD.y * avgE.y - 0.5 * E2)
             + (avgB.y * avgH.y - 0.5 * H2);
      ret.zz = (avgD.z * avgE.z - 0.5 * E2)
             + (avgB.z * avgH.z - 0.5 * H2);

      ret.xy = avgD.x*avgE.y + avgB.x*avgH.y;
      ret.yx = avgD.y*avgE.x + avgB.y*avgH.x;
      ret.xz = avgD.x*avgE.z + avgB.x*avgH.z;
      ret.zx = avgD.z*avgE.x + avgB.z*avgH.x;
      ret.yz = avgD.y*avgE.z + avgB.y*avgH.z;
      ret.zy = avgD.z*avgE.y + avgB.z*avgH.y;

      return ret;
    }
  };

public:

  /** Define AdvanceMethods using our stress tensor implementation. */
  template <typename Base>
  using AdvanceMethods = typename MBase::
      template AdvanceMethods<StressTensor, Base>;
};

} // namespace Integrate
} // namespace Output
} // namespace Fdtd

#endif // #ifndef FDTD_OUTPUT_INTEGRATE_MAXWELL_STRESS_TENSOR_HPP

