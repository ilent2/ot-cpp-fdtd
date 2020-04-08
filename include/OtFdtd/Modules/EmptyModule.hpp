/* Fdtd/EmptyModule.hpp - Base class for new decorator modules.
 *
 * Copyright (C) 2016 Isaac Lenton (aka ilent2)
 */

#ifndef FDTD_EMPTY_MODULE_HPP
#define FDTD_EMPTY_MODULE_HPP

namespace Fdtd {

/** Base class for Fdtd module groups. */
struct EmptyModule
{

  /** Describes the data stored at each field location. */
  template <typename Base> class ElementData : public Base {};

  /** Provides terminal cases for recursive methods.
   * Includes initialize and finalize methods.
   */
  template <typename Base> class Fallback : public Base {};

  /** Provides description of the simulation time steps.
   * Inherits from Fallback. */
  template <typename Base> class Timing : public Base {};

  /** Provides description of grid/dimensionality.
   * Inherits from Timing. */
  template <typename Base> class Indices : public Base {};

  /** Provides field data memory and initialization.
   * Inherits from Indices.  Elements allocated using ElementData. */
  template <typename Base> class FieldData : public Base {};

  /** Provides methods for each coordinate axis.
   * Inherits from FieldData. */
  template <typename Base> class CoordinateMethods : public Base {};

  /** Provides coordinate specific methods.
   * Inherits from CoordinateMethods. */
  template <typename Base> class Coordinates : public Base {};

  /** Provides methods for computing/storing field properties.
   * Inherits from Coordinates. */
  template <typename Base> class AdvanceMethods : public Base {};

  /** Provides methods for use in Evolve.
   * Inherits from AdvanceMethods. */
  template <typename Base> class EvolveMethods : public Base {};

  /** Implements methods for advancing/inspecting the simulation.
   * Inherits from  EvolveMethods. */
  template <typename Base> class Evolve : public Base {};

  /** Contains evolve methods.
   * Inherits from Evolve. */
  template <typename Base> class EvaluateMethods : public Base {};

  /** Contains the evaluate method called by the simulation object.
   * Inherits from EvaluateMethods. */
  template <typename Base> class Evaluate : public Base {};

};

} // namespace Fdtd

#endif // #ifndef FDTD_EMPTY_MODULE_HPP

