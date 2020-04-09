/* combined_index.hpp - Replicates the toolbox combined_index functions.
 *
 * Copyright (C) 2016 Isaac Lenton (aka ilent2)
 */

#ifndef TOOLBOX_COMBINED_INDEX_HPP
#define TOOLBOX_COMBINED_INDEX_HPP

#include <vector>
#include <tuple>

namespace Toolbox {

std::vector<unsigned> combined_index(
    const std::vector<int>& n, const std::vector<int>& m)
{
  std::vector<unsigned> ret(n.size());
  for (unsigned i = 0; i < n.size(); ++i)
  {
    ret[i] = n[i] * (n[i] + 1) + m[i];
  }

  return ret;
}

std::tuple<int, int> combined_index(unsigned ci)
{
  int n = floor(sqrt(ci));
  int m = ci - n*n - n;
  return std::make_tuple(n, m);
}

} // namespace Toolbox

#endif // #ifndef TOOLBOX_COMBINED_INDEX_HPP

