#include "rounding.hpp"

std::set<CGAL::Gmpq> build_nice_fraction_list(int max_denom) {
  // build a set of nice fractions in the interval (0,1],
  // like 1/10, 1/6, 1/3, 1/2, 2/3, etc.

  std::set<CGAL::Gmpq> set;
  set.insert(0);
  set.insert(1);
  for (auto denominator = 2; denominator <= max_denom; ++denominator)
    for (auto numerator = 1; numerator < denominator; ++numerator)
      set.insert(CGAL::Gmpq(numerator, denominator));
  return set;
}

auto find_nearest_nice_fraction(const CGAL::Gmpq &val,
                                const std::set<CGAL::Gmpq> &fracset) {

  // returns a pair: (nearest 'nice' fraction, distance)
  assert(fracset.size() >= 2);
  auto it = fracset.lower_bound(val);

  if (it == fracset.end())
    --it;

  auto dright = CGAL::abs(*it - val);

  if (it != fracset.begin()) {
    auto dleft = CGAL::abs(*(--it) - val);
    if (dleft < dright)
      return std::make_pair(*it, dleft);
    else
      it++;
  }
  return std::make_pair(*it, dright);
}

std::pair<CGAL::Gmpq, bool>
roundNice(const CGAL::Gmpq& x, const std::set<CGAL::Gmpq>& fracset, int digit_threshold, double eps) {

  // round the input to a "nice" fraction if it is eps-close (in absolute sense)
  // to such a fraction

  auto n = std::max(x.numerator().approximate_decimal_length(),
                    x.denominator().approximate_decimal_length());
  if (n < digit_threshold)
    return {x, false};

  CGAL::Gmpz integral_part(std::floor(x.to_double()));
  auto frac = x - integral_part;
  auto res = find_nearest_nice_fraction(CGAL::abs(frac),fracset);

  if (res.second.to_double() < eps) {
    return {integral_part + CGAL::sign(frac) * res.first, true};
  }
  return {x, false};
}

std::pair<Polygon_2, bool> round_vertices(const Polygon_2 &P,
                                          const std::set<CGAL::Gmpq> &fracset,
                                          int digit_threshold, double error) {
  Polygon_2 res;
  bool did_we_round = false;
  for (auto i = P.vertices_begin(); i != P.vertices_end(); ++i) {
    auto p = *i;
    auto X = roundNice(p.x(), fracset, digit_threshold, error);
    auto Y = roundNice(p.y(), fracset, digit_threshold, error);
    res.push_back(Point_2(X.first, Y.first));
    did_we_round |= (X.second || Y.second);
  }
  return {res, did_we_round};
}
