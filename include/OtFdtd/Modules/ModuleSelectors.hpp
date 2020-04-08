/* ModuleSelectors.hpp - Selectors for use with DecoratedChild.
 *
 * Copyright (C) 2016 Isaac Lenton (aka ilent2)
 */

#ifndef FDTD_MODULE_SELECTORS_HPP
#define FDTD_MODULE_SELECTORS_HPP

namespace Fdtd {

/** Selector for element data groups. */
template <typename Parent> struct ElementDataSelector
{ template <typename Base> using Type
    = typename Parent::template ElementData<Base>; };

/** Selector for fallbacks: terminals cases for traversal functions. */
template <typename Parent> struct FallbackSelector
{ template <typename Base> using Type
    = typename Parent::template Fallback<Base>; };

/** Selector for Timing data. */
template <typename Parent> struct TimingSelector
{ template <typename Base> using Type
    = typename Parent::template Timing<Base>; };

/** Selector for Indices data. */
template <typename Parent> struct IndicesSelector
{ template <typename Base> using Type
    = typename Parent::template Indices<Base>; };

/** Selector for FieldData. */
template <typename Parent> struct FieldDataSelector
{ template <typename Base> using Type
    = typename Parent::template FieldData<Base>; };

/** Selector for CoordinateMethods. */
template <typename Parent> struct CoordinateMethodsSelector
{ template <typename Base> using Type
    = typename Parent::template CoordinateMethods<Base>; };

/** Selector for Coordinates. */
template <typename Parent> struct CoordinatesSelector
{ template <typename Base> using Type
    = typename Parent::template Coordinates<Base>; };

/** Selector for AdvanceMethods. */
template <typename Parent> struct AdvanceMethodsSelector
{ template <typename Base> using Type
    = typename Parent::template AdvanceMethods<Base>; };

/** Selector for EvolveMethods. */
template <typename Parent> struct EvolveMethodsSelector
{ template <typename Base> using Type
    = typename Parent::template EvolveMethods<Base>; };

/** Selector for Evolve. */
template <typename Parent> struct EvolveSelector
{ template <typename Base> using Type
    = typename Parent::template Evolve<Base>; };

/** Selector for EvaluateMethods. */
template <typename Parent> struct EvaluateMethodsSelector
{ template <typename Base> using Type
    = typename Parent::template EvaluateMethods<Base>; };

/** Selector for Evaluate. */
template <typename Parent> struct EvaluateSelector
{ template <typename Base> using Type
    = typename Parent::template Evaluate<Base>; };

} // namespace Fdtd

#endif // #ifndef FDTD_MODULE_SELECTORS_HPP

