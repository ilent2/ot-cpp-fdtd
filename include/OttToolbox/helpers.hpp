/* helpers.hpp - All the other functions we could want.
 *
 * Copyright (C) 2016 Isaac Lenton (aka ilent2)
 */

#ifndef TOOLBOX_HELPERS_HPP
#define TOOLBOX_HELPERS_HPP

namespace Toolbox {

/** Replace NaNs with zero for Re(z) == NaN || Im(z) == NaN. */
void StripNans(Utilities::Vec3c& a)
{
  if (std::isnan(std::real(a.x)) || std::isnan(std::imag(a.x))) a.x = 0.0;
  if (std::isnan(std::real(a.y)) || std::isnan(std::imag(a.y))) a.y = 0.0;
  if (std::isnan(std::real(a.z)) || std::isnan(std::imag(a.z))) a.z = 0.0;
}

/** Reproduce the toolbox function rtpv2xyzv. */
Utilities::Vec3c rtpv2xyzv(const Utilities::VecRtp& rtp,
    const Utilities::VecRtpc& inp)
{
  const auto& rv = inp.r;
  const auto& tv = inp.t;
  const auto& pv = inp.p;

  double theta = rtp.t;
  double phi = rtp.p;

  Utilities::Vec3c ret;
  ret.x = rv*sin(theta)*cos(phi) + tv*cos(theta)*cos(phi) - pv*sin(phi);
  ret.y = rv*sin(theta)*sin(phi) + tv*cos(theta)*sin(phi) + pv*cos(phi);
  ret.z = rv*cos(theta) - tv*sin(theta);

  return ret;
}

} // namespace Toolbox

#endif // #ifndef TOOLBOX_HELPERS_HPP

