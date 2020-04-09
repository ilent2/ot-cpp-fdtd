/* *********************** BEGIN LICENSE BLOCK *******************************\
 *
 * Implementation of Optics Toolbox electromagnetic_field_xyz function.
 * Copyright (C) 2016 Isaac Lenton (aka ilent2)
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 *
\* *********************** END LICENSE BLOCK *********************************/

// This version is not optimised, it simply tries to implement the desired
// functionality by following the original optical toolbox implementation.
//
// We've tried to clean it up a bit since the first version, but
// the above clause still applies, it's still not optimised.

#ifndef TOOLBOX_EM_FIELD_XYZ_HPP
#define TOOLBOX_EM_FIELD_XYZ_HPP

#include <tuple>
#include <vector>
#include <complex>

// Calculation of the E/H fields requires the VSWF, which in turn requires
// calculation of the Legendre function and spherical bessil functions.
// We use the already implemented C++ Optical Toolbox classes for these.
#include "bessel.hpp"
#include "vsh.hpp"

// We use the Utilities types for the FDTD simulation instead of our own
#include "Utilities/Vec.hpp"
#include "Utilities/VecRtp.hpp"

// Include the other files with our Toolbox namespace
#include "helpers.hpp"
#include "vswf.hpp"

namespace Toolbox {

/** Calculate the electromagnetic field values using VSWF.
 *
 * We've split the function into two since the last revision,
 * this version iterates over locations rather than unique n values.
 *
 * @param locs  : Locations to calculate the field at.
 * @param m     : A vector or +/- 1 describing the component contributions.
 * @param n     : A vector of positive integers for component indices.
 * @param a     : The complex amplitudes.
 * @param b     : The complex amplitudes.
 * @param wavenumber  : The wavenumber in the background medium.
 *
 * m, n, a and b should be the same length.
 *
 * @return A vector of complex E field values.
 */
std::vector<Utilities::Vec3c> emFieldXyzE(
    const std::vector<Utilities::Vec3d>& locs,
    const std::vector<int>& m, const std::vector<int>& n,
    const std::vector<std::complex<double>>& a,
    const std::vector<std::complex<double>>& b,
    double wavenumber)
{
  std::vector<Utilities::Vec3c> ret;
  ret.reserve(locs.size());

  // What's out maximum n
  int nmax = *std::max_element(n.begin(), n.end());

  // Prepare to calculate the Spherical Bessel functions for the VSWFs
  ott::bessel<double> bessel(nmax);

  // Prepare to calculate spherical harmonics
  Toolbox::Vsh<double> vsh(nmax);

  for (const auto& loc : locs)
  {
    // Convert the coordinate to a spherical coordinate
    Utilities::VecRtp rtp(loc);

    // Declare the memory for storing the result (in spherical coordinates)
    Utilities::VecRtpc value;

    // Update bessel/vsh with our coordinate
    bessel.change_x(rtp.r * wavenumber);
    vsh.change_x(rtp.t, rtp.p);

    for (unsigned i = 0; i < m.size(); ++i)
    {
      Utilities::VecRtpc M3, N3;
      std::tie(M3, N3) = vswf(n[i], m[i], wavenumber, rtp, vsh, bessel);
      value = value.add_components(M3.mul_components(a[i]));
      value = value.add_components(N3.mul_components(b[i]));
    }

    // Convert result to Cartesian coordinate
    ret.push_back(rtpv2xyzv(rtp, value));

    // Remove NaNs
    StripNans(ret.back());
  }

  return ret;
}

/** Calculate the electromagnetic field values using VSWF.
 *
 * We've split the function into two since the last revision,
 * this version iterates over locations rather than unique n values.
 *
 * @param locs  : Locations to calculate the field at.
 * @param m     : A vector or +/- 1 describing the component contributions.
 * @param n     : A vector of positive integers for component indices.
 * @param a     : The complex amplitudes.
 * @param b     : The complex amplitudes.
 * @param wavenumber  : The wavenumber in the background medium.
 *
 * m, n, a and b should be the same length.
 *
 * @return A vector of complex H field values.
 *    The H field values have not been multiplied by the impedance.
 */
std::vector<Utilities::Vec3c> emFieldXyzH(
    const std::vector<Utilities::Vec3d>& locs,
    const std::vector<int>& m, const std::vector<int>& n,
    const std::vector<std::complex<double>>& a,
    const std::vector<std::complex<double>>& b,
    double wavenumber)
{
  std::vector<Utilities::Vec3c> ret;
  ret.reserve(locs.size());

  // What's out maximum n
  int nmax = *std::max_element(n.begin(), n.end());

  // Prepare to calculate the Spherical Bessel functions for the VSWFs
  ott::bessel<double> bessel(nmax);

  // Prepare to calculate spherical harmonics
  Toolbox::Vsh<double> vsh(nmax);

  for (const auto& loc : locs)
  {
    // Convert the coordinate to a spherical coordinate
    Utilities::VecRtp rtp(loc);

    // Declare the memory for storing the result (in spherical coordinates)
    Utilities::VecRtpc value;

    // Update bessel/vsh with our coordinate
    bessel.change_x(rtp.r * wavenumber);
    vsh.change_x(rtp.t, rtp.p);

    for (unsigned i = 0; i < m.size(); ++i)
    {
      Utilities::VecRtpc M3, N3;
      std::tie(M3, N3) = vswf(n[i], m[i], wavenumber, rtp, vsh, bessel);
      value = value.add_components(M3.mul_components(b[i]));
      value = value.add_components(N3.mul_components(a[i]));
    }

    value = value.mul_components(std::complex<double>(0.0, 1.0));

    // Convert result to Cartesian coordinate
    ret.push_back(rtpv2xyzv(rtp, value));

    // Remove NaNs
    StripNans(ret.back());
  }

  return ret;
}

/** Calculate the electromagnetic field values using VSWF.
 *
 * We've split the function into two since the last revision,
 * this version iterates over locations rather than unique n values.
 *
 * @param locs  : Locations to calculate the field at.
 * @param m     : A vector or +/- 1 describing the component contributions.
 * @param n     : A vector of positive integers for component indices.
 * @param a     : The complex amplitudes.
 * @param b     : The complex amplitudes.
 * @param wavenumber  : The wavenumber in the background medium.
 *
 * TODO: This shares a lot in common with the other emFieldXyz
 *    methods, perhaps we should merge them again?
 *
 * m, n, a and b should be the same length.
 *
 * @return A tuple with vectors of complex E and H field values.
 *    The H field values have not been multiplied by the impedance.
 */
std::tuple<std::vector<Utilities::Vec3c>, std::vector<Utilities::Vec3c>>
emFieldXyz(
    const std::vector<Utilities::Vec3d>& locs,
    const std::vector<int>& m, const std::vector<int>& n,
    const std::vector<std::complex<double>>& a,
    const std::vector<std::complex<double>>& b,
    double wavenumber)
{
  std::vector<Utilities::Vec3c> retE;
  std::vector<Utilities::Vec3c> retH;
  retE.reserve(locs.size());
  retH.reserve(locs.size());

  // What's out maximum n
  int nmax = *std::max_element(n.begin(), n.end());

  // Prepare to calculate the Spherical Bessel functions for the VSWFs
  ott::bessel<double> bessel(nmax);

  // Prepare to calculate spherical harmonics
  Toolbox::Vsh<double> vsh(nmax);

  for (const auto& loc : locs)
  {
    // Convert the coordinate to a spherical coordinate
    Utilities::VecRtp rtp(loc);

    // Declare the memory for storing the result (in spherical coordinates)
    Utilities::VecRtpc valueE;
    Utilities::VecRtpc valueH;

    // Update bessel/vsh with our coordinate
    bessel.change_x(rtp.r * wavenumber);
    vsh.change_x(rtp.t, rtp.p);

    for (unsigned i = 0; i < m.size(); ++i)
    {
      Utilities::VecRtpc M3, N3;
      std::tie(M3, N3) = vswf(n[i], m[i], wavenumber, rtp, vsh, bessel);

      valueE = valueE.add_components(M3.mul_components(a[i]));
      valueE = valueE.add_components(N3.mul_components(b[i]));

      valueH = valueH.add_components(M3.mul_components(b[i]));
      valueH = valueH.add_components(N3.mul_components(a[i]));
    }

    valueH = valueH.mul_components(std::complex<double>(0.0, 1.0));

    // Convert result to Cartesian coordinate
    retE.push_back(rtpv2xyzv(rtp, valueE));
    retH.push_back(rtpv2xyzv(rtp, valueH));

    // Remove NaNs
    StripNans(retE.back());
    StripNans(retH.back());
  }

  return std::make_tuple(retE, retH);
}

} /* namespace Toolbox */

#endif // #ifndef TOOLBOX_EM_FIELD_XYZ_HPP

