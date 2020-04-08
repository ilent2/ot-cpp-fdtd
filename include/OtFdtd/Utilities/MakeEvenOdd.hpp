/** Helper functions for making things even and odd
 *
 *
 * Part of ot-cpp-fdtd, Copyright 2020 Isaac Lenton
 * See LICENSE for details about using/distributing this file.
 */

constexpr unsigned make_odd(unsigned i) { return 2*(i/2) + 1; }
constexpr unsigned make_even(unsigned i) { return 2*(i/2); }

