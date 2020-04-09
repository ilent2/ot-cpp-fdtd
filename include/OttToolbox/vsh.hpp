/* Toolbox/vsh.hpp - An interface to Legendre that calculates VSH.
 *
 * This is similar to the toolbox function vsh.m
 *
 * Copyright (C) 2016 Isaac Lenton (aka ilent2)
 */

#ifndef TOOLBOX_VSH_HPP
#define TOOLBOX_VSH_HPP

#include "legendre.hpp"
#include <cmath>
#include <complex>
#include <tuple>
#include "Utilities/VecRtp.hpp"

namespace Toolbox {

/** A wrapper for legendre providing the vector spherical harmonics.
 *
 * TODO: Generalize types
 */
template <typename T>
class Vsh : public ott::legendre<double>
{

  static_assert(std::is_same<T, double>::value,
      "Only implemented for double.");

  double phi_now;

public:

  /** Construct a legendre instance. */
  Vsh(int nmax) : legendre(nmax, 0.0) {}

  /** Provide a new change location method with theta and phi. */
  void change_x(double theta, double phi)
  {
    legendre::change_x(theta);
    phi_now = phi;
  }

  /** Produces similar output to the toolbox spharm.m function. */
  std::tuple<std::complex<double>, std::complex<double>, std::complex<double>>
  spharm(int n, int m)
  {
    double lY = legendre::Y(n, m);
    double ltY = legendre::dtY(n, m);
    double lpY = legendre::dpY(n, m);

    std::complex<double> expphi = std::exp(
        std::complex<double>(0.0, m) * phi_now);

    return std::make_tuple(lY * expphi, ltY * expphi,
        std::complex<double>(0.0, 1.0) * lpY * expphi);
  }

  /** Produces similar output to the toolbox vsh.m function. */
  std::tuple<Utilities::VecRtpc, Utilities::VecRtpc, Utilities::VecRtpc>
  vsh(int n, int m)
  {
    std::complex<double> Y, dtY, dpY;
    std::tie(Y, dtY, dpY) = spharm(n, m);

    return std::make_tuple(Utilities::VecRtpc(0.0, dtY, dpY),
                           Utilities::VecRtpc(0.0, dpY, -dtY),
                           Utilities::VecRtpc(Y, 0.0, 0.0));
  }
};

} // namespace Toolbox

#endif // #ifndef TOOLBOX_VSH_HPP

