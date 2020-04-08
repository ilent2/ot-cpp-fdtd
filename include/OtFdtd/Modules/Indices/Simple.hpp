/* Fdtd/Indices/Simple.hpp - Simple indices for a rect-linear grid.
 *
 * TODO: Generalize for 1-D and 2-D again.
 *
 * Copyright (C) 2016 Isaac Lenton (aka ilent2)
 */

#ifndef FDTD_INDICES_SIMPLE_HPP
#define FDTD_INDICES_SIMPLE_HPP

#include "Utilities/DecoratedChild.hpp"
#include "../ModuleSelectors.hpp"
#include "../EmptyModule.hpp"
#include "Pbc.hpp"
#include "Linear.hpp"

namespace Fdtd {
namespace Indices {


template <unsigned x, unsigned y, unsigned z>
struct Simple : public EmptyModule
{

  template <typename Base>
  using Fallback = typename Utilities::DecoratedChild<FallbackSelector,
      Linear<0, z>, Linear<1, y>, Linear<2, x>>::template Type<Base>;

  template <typename Base>
  using Indices = typename Utilities::DecoratedChild<IndicesSelector,
      Linear<0, z>, Linear<1, y>, Linear<2, x>>::template Type<Base>;

};

/** Determine if base is a 1-D Cartesian coordinate system (i.e no PBC). */
template <typename Base>
struct is_cartesian_1d
{
  static constexpr bool value = Base::AxisIndex == 1
      && !has_pbc<Base, 0>::value;
};

/** Determine if base is a 2-D Cartesian coordinate system (i.e no PBC). */
template <typename Base>
struct is_cartesian_2d
{
  static constexpr bool value = Base::AxisIndex == 2
      && !has_pbc<Base, 0>::value && !has_pbc<Base, 1>::value;
};

/** Determine if base is a 3-D Cartesian coordinate system (i.e no PBC). */
template <typename Base>
struct is_cartesian_3d
{
  static constexpr bool value = Base::AxisIndex == 3
      && !has_pbc<Base, 0>::value
      && !has_pbc<Base, 1>::value
      && !has_pbc<Base, 2>::value;
};


} // namespace Indices
} // namespace Fdtd

#endif // #ifndef FDTD_INDICES_SIMPLE_HPP

