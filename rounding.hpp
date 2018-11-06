#ifndef ROUNDINGHPP
#define ROUNDINGHPP

#include <utility>
#include <set>
#include "polygon.hpp"

std::set<CGAL::Gmpq> build_nice_fraction_list(int max_denom);
std::pair<Polygon_2, bool> round_vertices(const Polygon_2 &P,
                                          const std::set<CGAL::Gmpq> &fracset,
                                          int digit_threshold, double error);

#endif


