/* Utilities/DecoratedLayers.hpp - Decorate object with layers of children.
 *
 * Uses the Utilities/DecoratedChild class for each child selector.
 *
 * Copyright (C) 2016 Isaac Lenton (aka ilent2)
 */

#ifndef UTILITIES_DECORATED_LAYERS_HPP
#define UTILITIES_DECORATED_LAYERS_HPP

#include "DecoratedChild.hpp"

namespace Utilities {

/** Creates a new template type taking a single template base class.
 * Performs a similar operation to @see DecoratedChild but repeats
 * the operation for each child selector.
 *
 * With only 1 child, this class is equivalent to DecoratedChild.
 *
 * @tparam Children : Template types describing which child member to use.
 * @tparam Decorators : Classes providing the child members.
 *
 * Usage:
 *
 * template <typename Child> struct Sel1
 * { template <typename Base> using Type = Child::TheChildType1<Base>; };
 *
 * template <typename Child> struct Sel1
 * { template <typename Base> using Type = Child::TheChildType2<Base>; };
 *
 * struct Dec1 {
 *   template <typename Base> class TheChildType1 {};
 *   template <typename Base> class TheChildType2 {};
 * };
 *
 * struct Dec2 {
 *   template <typename Base> class TheChildType1 {};
 *   template <typename Base> class TheChildType2 {};
 * };
 *
 * template <typename Base>
 * using NewType = DecoratedLayers<Sel1, Sel2, Dec1, Dec2>::Type<Base>;
 *
 * Would be equivalent to:
 *
 * template <typename Base>
 * using NewType = Dec2::TheChildType2<Dec1::TheChildType2<
 *                 Dec2::TheChildType1<Dec1::TheChildType1<Base>>>>;
 */
template <template <typename> class... Children>
class DecoratedLayers
{
public:

  template <typename... Decorators>
  class For
  {

    template <template <typename> class First,
        template <typename> class... Others>
    struct Assemble
    {
      template <typename Base>
      using Type = typename Assemble<Others...>::template Type<
          typename DecoratedChild<First, Decorators...>::template Type<Base>>;
    };

    template <template <typename> class Last>
    struct Assemble<Last>
    {
      template <typename Base>
      using Type = typename DecoratedChild<Last,
          Decorators...>::template Type<Base>;
    };

  public:

    /** The resulting type. */
    template <typename Base>
    using Type = typename Assemble<Children...>::template Type<Base>;

  };

};

} // namepsace Utilities

#endif // #ifndef UTILITIES_DECORATED_LAYERS_HPP

