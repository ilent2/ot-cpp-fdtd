/* Utilities/MinMaxAttribute.hpp - Calculate the min/max attribute value.
 *
 * Copyright (C) 2016 Isaac Lenton (aka ilent2)
 */

#ifndef UTILITIES_MIN_MAX_ATTRIBUTE_HPP
#define UTILITIES_MIN_MAX_ATTRIBUTE_HPP

namespace Utilities {

/** Calculate the minimum and maximum wavelengths of supplied materials. */
template <template <typename> class Attr, typename First, typename... Types>
struct MinMaxAttribute
{
  static constexpr double min = std::min(Attr<First>::value,
      MinMaxAttribute<Attr, Types...>::min);
  static constexpr double max = std::max(Attr<First>::value,
      MinMaxAttribute<Attr, Types...>::max);
};

/** Terminal case for recursive MinMaxWavelength. */
template <template <typename> class Attr, typename Last>
struct MinMaxAttribute<Attr, Last>
{
  static constexpr double min = Attr<Last>::value;
  static constexpr double max = Attr<Last>::value;
};

} // namespace Utilities

#endif // #ifndef UTILITIES_MIN_MAX_ATTRIBUTE_HPP

