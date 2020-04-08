/* Utilities/DecoratedChild.hpp - Decorate a object with object child.
 *
 * Copyright (C) 2016 Isaac Lenton (aka ilent2)
 */

#ifndef UTILITIES_DECORATED_CHILD_HPP
#define UTILITIES_DECORATED_CHILD_HPP

namespace Utilities {

/** Creates a new template type taking a single template base class.
 * The base is decorated with the member classes, specified by the
 * Child template parameter, from the decorator class list.
 *
 * @tparam Child : Describes the child class to use from Decorators.
 * @tparam Decorators : The objects containing the decorating classes.
 *
 * Usage:
 *
 * template <typename Child> struct ChildSelector
 * { template <typename Base> using Type = Child::TheChildType<Base>; };
 *
 * struct Dec1 { template <typename Base> class TheChildType {}; };
 * struct Dec2 { template <typename Base> class TheChildType {}; };
 * struct Dec3 { template <typename Base> class TheChildType {}; };
 *
 * template <typename Base>
 * using NewType = DecoratedChild<ChildSelector, Dec1, Dec2, Dec3>::Type<Base>;
 *
 * Would be equivalent to:
 *
 * template <typename Base>
 * using NewType = Dec3::TheChildType<
 *    Dec2::TheChildType<Dec1::TheChildType<Base>>>;
 */
template <template <typename> class Child, typename... Decorators>
class DecoratedChild
{

  template <typename First, typename... Others>
  struct Assemble
  {
    template <typename Base>
    using Type = typename Assemble<Others...>::template Type<
        typename Child<First>::template Type<Base>>;
  };

  template <typename Last>
  struct Assemble<Last>
  {
    template <typename Base>
    using Type = typename Child<Last>::template Type<Base>;
  };

public:

  /** The resulting type. */
  template <typename Base>
  using Type = typename Assemble<Decorators...>::template Type<Base>;

};

} // namespace Utilities

#endif // #ifndef UTILITIES_DECORATED_CHILD_HPP

