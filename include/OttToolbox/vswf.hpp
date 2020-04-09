/* Toolbox/vswf.hpp - Implementation of vswf.m
 *
 * Copyright (C) 2016 Isaac Lenton (aka ilent2)
 */

#ifndef TOOLBOX_VSWF_HPP
#define TOOLBOX_VSWF_HPP

#include <tuple>
#include "Utilities/VecRtp.hpp"
#include "vsh.hpp"

namespace Toolbox {

std::tuple<Utilities::VecRtpc, Utilities::VecRtpc> vswf(
    int n, int m, double wavenumber, const Utilities::VecRtp& rtp,
    Vsh<double>& vsh, const ott::bessel<double>& bessel)
{
  Utilities::VecRtpc M, N;

  if (n == 0) return std::make_tuple(M, N);
  double Nn = pow(n*(n+1.0), -0.5);

  Utilities::VecRtpc B, C, P;
  std::tie(B, C, P) = vsh.vsh(n, m);

  double besseljn = bessel.j()[n];
  double besseljn1 = bessel.j()[n-1];
  double kr = rtp.r * wavenumber;

  M = C.mul_components(Nn * besseljn);

  // TODO: Hmm, perhaps a tolerance should be used here?
  if (kr == 0.0)
  {
    N = n != 1 ? Utilities::VecRtpc()
        : (P.add_components(B)).mul_components((2.0/3.0) * Nn);
  }
  else
  {
    N = P.mul_components(Nn * (n*(n+1.0)/kr * besseljn)).add_components(
        B.mul_components(Nn * (besseljn1 - besseljn*n/kr)));
  }

  return std::make_tuple(M, N);
}

} // namespace Toolbox

#endif // #ifndef TOOLBOX_VSWF_HPP

