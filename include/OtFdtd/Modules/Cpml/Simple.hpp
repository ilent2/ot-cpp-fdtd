/* Fdtd/Cpml/Simple.hpp - Automatically deduced CPML.
 *
 * A CPML type suitable for 1, 2 and 3 dimensional Cartesian,
 * polar, cylindrical and spherical coordinates.
 *
 * The type of CPML to use is determined using the Indices
 * methods for querying type in Indices/Polar and Indices/Simple.
 *
 * TODO: Generalize for variable time steps.
 *
 * Copyright (C) 2016 Isaac Lenton (aka ilent2)
 */

#ifndef FDTD_CPML_SIMPLE_HPP
#define FDTD_CPML_SIMPLE_HPP

#include "DefaultTypes.hpp"             // default parameters for Simple
#include "../EmptyModule.hpp"           // base class for Simple
#include "../Indices/Polar.hpp"         // for coordinate type query
#include "../Indices/Simple.hpp"        // for coordinate type query
#include "../Coordinate/Stretched.hpp"  // for stretched coordinates
#include "../ModuleSelectors.hpp"       // for Assemble class
#include "CoordinateFromStretching.hpp" // for stretched coordinates
#include "Cpml.hpp"                     // evolve methods definition

namespace Fdtd {
namespace Cpml {

/** Internal methods used by Fdtd::Cpml::Simple. */
namespace _Simple {

/** Implementation of CPML with stretched coordinates for one end of a axis.
 *
 * @tparam Depth        Depth of CPML layer from simulation boundary.
 * @tparam Axis         Axis to add CPML to.
 * @tparam Forward      End to apply CPML to (true: start, false: end).
 * @tparam eps          Permittivity inside CPML.
 * @tparam mu           Permeability inside CPML.
 * @tparam Grading      Grading used for CPML parameters.
 * @tparam Stretching   Stretching used for stretched coordinates.
 */
template <unsigned Depth, int tAxis, bool Forward,
    const double& TimeStep, const double& eps, const double& mu,
    typename Grading = DefaultStretching,
    typename Stretching = DefaultStretching>
struct Axis
{
  template <typename Base>
  using CoordinateMethods = typename Coordinate::Stretched<tAxis,
      CoordinateFromStretching<Stretching, Depth, Forward>
      >::template CoordinateMethods<Base>;

  template <typename Base>
  using EvolveMethods = Cpml<Base, Depth, tAxis, Forward, TimeStep,
      eps, mu, Grading>;
};

/** Assemble the Coordinate and EvolveMethods members. */
template <typename... Parts>
struct Assemble
{
  template <typename Base>
  using CoordinateMethods = typename Utilities::DecoratedChild<
      CoordinateMethodsSelector, Parts...>::template Type<Base>;

  template <typename Base>
  using EvolveMethods = typename Utilities::DecoratedChild<
      EvolveMethodsSelector, Parts...>::template Type<Base>;
};

/** Implementation of CPML with stretched coordinates for a axis. */
template <unsigned Depth, int tAxis,
    const double& TimeStep, const double& eps, const double& mu,
    typename Grading = DefaultStretching,
    typename Stretching = DefaultStretching>
using FullAxis = Assemble<
    Axis<Depth, tAxis, true, TimeStep, eps, mu, Grading, Stretching>,
    Axis<Depth, tAxis, false, TimeStep, eps, mu, Grading, Stretching>>;

/** Signature for selecting appropriate Coordinate/EvolveMethods. */
template <typename Base, unsigned Depth,
    const double& TimeStep, const double& eps, const double& mu,
    typename Grading, typename Stretching, typename Enable = void>
class Methods
{
  // Coordinate/EvolveMethods unimplemented, prevent direct instantiation
};

/** 1-D Cartesian */
template <typename Base, unsigned Depth,
    const double& TimeStep, const double& eps, const double& mu,
    typename Grading, typename Stretching>
class Methods<Base, Depth, TimeStep, eps, mu, Grading, Stretching,
    std::enable_if_t<Indices::is_cartesian_1d<Base>::value>>
: public FullAxis<Depth, 0, TimeStep, eps, mu, Grading, Stretching>
{
};

/** 2-D Cartesian */
template <typename Base, unsigned Depth,
    const double& TimeStep, const double& eps, const double& mu,
    typename Grading, typename Stretching>
class Methods<Base, Depth, TimeStep, eps, mu, Grading, Stretching,
    std::enable_if_t<Indices::is_cartesian_2d<Base>::value>>
: public Assemble<
    FullAxis<Depth, 0, TimeStep, eps, mu, Grading, Stretching>,
    FullAxis<Depth, 1, TimeStep, eps, mu, Grading, Stretching>>
{
};

/** 3-D Cartesian */
template <typename Base, unsigned Depth,
    const double& TimeStep, const double& eps, const double& mu,
    typename Grading, typename Stretching>
class Methods<Base, Depth, TimeStep, eps, mu, Grading, Stretching,
    std::enable_if_t<Indices::is_cartesian_3d<Base>::value>>
: public Assemble<
    FullAxis<Depth, 0, TimeStep, eps, mu, Grading, Stretching>,
    FullAxis<Depth, 1, TimeStep, eps, mu, Grading, Stretching>,
    FullAxis<Depth, 2, TimeStep, eps, mu, Grading, Stretching>>
{
};

/** 2-D Polar */
template <typename Base, unsigned Depth,
    const double& TimeStep, const double& eps, const double& mu,
    typename Grading, typename Stretching>
class Methods<Base, Depth, TimeStep, eps, mu, Grading, Stretching,
    std::enable_if_t<Indices::is_polar<Base>::value>>
: public Axis<Depth, 1, true, TimeStep, eps, mu, Grading, Stretching>
{
};

/** 3-D Cylindrical */
template <typename Base, unsigned Depth,
    const double& TimeStep, const double& eps, const double& mu,
    typename Grading, typename Stretching>
class Methods<Base, Depth, TimeStep, eps, mu, Grading, Stretching,
    std::enable_if_t<Indices::is_cylindrical<Base>::value>>
: public Assemble<
    Axis<Depth, 1, true, TimeStep, eps, mu, Grading, Stretching>,
    FullAxis<Depth, 2, TimeStep, eps, mu, Grading, Stretching>>
{
};

/** 3-D Spherical */
template <typename Base, unsigned Depth,
    const double& TimeStep, const double& eps, const double& mu,
    typename Grading, typename Stretching>
class Methods<Base, Depth, TimeStep, eps, mu, Grading, Stretching,
    std::enable_if_t<Indices::is_spherical<Base>::value>>
: public Axis<Depth, 2, true, TimeStep, eps, mu, Grading, Stretching>
{
};

} // namespace _Simple

/** Automatically deduced CPML.
 *
 * A CPML type suitable for 1, 2 and 3 dimensional Cartesian,
 * polar, cylindrical and spherical coordinates.
 *
 * This used to be the signature class, with an Enable template
 * parameter, however our new model moves the Base template parameter
 * inside the module group classes, so we need to move Enable too.
 *
 * @tparam Depth        Depth of CPML layer from simulation boundary.
 * @tparam TimeStep     A parameter used by the CPML.
 * @tparam Grading      Grading used for CPML parameters.
 * @tparam Stretching   Stretching used for stretched coordinates.
 */
template <unsigned Depth,
    const double& TimeStep, const double& eps, const double& mu,
    typename Grading = DefaultStretching,
    typename Stretching = DefaultStretching>
class Simple : public EmptyModule
{
  /** Select the appropriate CPML implementation for the given Base. */
  template <typename Base>
  using Modules = typename _Simple::Methods<Base,
      Depth, TimeStep, eps, mu, Grading, Stretching>;

public:

  template <typename Base> using CoordinateMethods =
      typename Modules<Base>::template CoordinateMethods<Base>;

  template <typename Base>
  using EvolveMethods = typename Modules<Base>::template EvolveMethods<Base>;
};

/** Apply CPML to a single axis.
 *
 * A CPML type suitable for 1, 2 and 3 dimensional Cartesian coordinates
 * and possible some other configurations.
 *
 * @tparam Depth        Depth of CPML layer from simulation boundary.
 * @tparam Axis         The axis to apply the CPML to.
 * @tparam TimeStep     A parameter used by the CPML.
 * @tparam Permittivity Vacuum permittivity.
 * @tparam Permeability Vacuum permeability.
 * @tparam Grading      Grading used for CPML parameters.
 * @tparam Stretching   Stretching used for stretched coordinates.
 */
template <unsigned Depth, int Axis,
    const double& TimeStep, const double& eps, const double& mu,
    typename Grading = DefaultStretching,
    typename Stretching = DefaultStretching>
class FullAxis : public EmptyModule
{
  /** Select the appropriate CPML implementation for the given Base. */
  template <typename Base>
  using Modules = typename _Simple::FullAxis<
      Depth, Axis, TimeStep, eps, mu, Grading, Stretching>;

public:

  template <typename Base> using CoordinateMethods =
      typename Modules<Base>::template CoordinateMethods<Base>;

  template <typename Base>
  using EvolveMethods = typename Modules<Base>::template EvolveMethods<Base>;
};

} // namespace Cpml
} // namespace Fdtd

#endif // #ifndef FDTD_CPML_SIMPLE_HPP

