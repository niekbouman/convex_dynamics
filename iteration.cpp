// standard includes
#include <iostream>
#include <utility>

#include <fstream>
#include <cassert>
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <set>
#include <utility>
// includes for defining the Voronoi diagram adaptor

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Boolean_set_operations_2.h>

//#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Voronoi_diagram_2.h>
#include <CGAL/Delaunay_triangulation_adaptation_traits_2.h>
#include <CGAL/Delaunay_triangulation_adaptation_policies_2.h>
//
#include <CGAL/Polygon_2.h>
//#include <CGAL/basic.h>
#include <CGAL/Minkowski_sum_2.h>
#include <CGAL/Direction_2.h>

//#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/ch_graham_andrew.h>
#include <CGAL/ch_melkman.h>
#include <CGAL/Aff_transformation_2.h>

#include "json.hpp" // nlohmann's json lib

// typedefs for defining the adaptor
// typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
// #include <CGAL/Simple_homogeneous.h>
#include <CGAL/Gmpz.h>
#include <CGAL/Gmpq.h>

typedef CGAL::Simple_cartesian<CGAL::Gmpq> K; 

//typedef CGAL::Exact_predicates_exact_constructions_kernel K;
typedef CGAL::Delaunay_triangulation_2<K> DT;
typedef CGAL::Delaunay_triangulation_adaptation_traits_2<DT> AT;
typedef CGAL::Delaunay_triangulation_caching_degeneracy_removal_policy_2<DT> AP;
typedef CGAL::Voronoi_diagram_2<DT, AT, AP> VD;

// typedef for the result type of the point location
typedef AT::Site_2 Site_2;
typedef AT::Point_2 Point_2;
typedef VD::Locate_result Locate_result;
typedef VD::Vertex_handle Vertex_handle;
typedef VD::Face_handle Face_handle;
typedef VD::Halfedge_handle Halfedge_handle;
typedef VD::Ccb_halfedge_circulator Ccb_halfedge_circulator;

typedef CGAL::Aff_transformation_2<K> Transformation;
typedef CGAL::Polygon_2<K> Polygon_2;
typedef CGAL::Vector_2<K> Vector_2;
typedef CGAL::Direction_2<K> Direction_2;
typedef CGAL::Polygon_with_holes_2<K> Polygon_with_holes_2;
typedef std::list<Polygon_with_holes_2> Pwh_list_2;

auto build_nice_fraction_list()
{
  // build a set of nice fractions in the interval (0,1],
  // like 1/10, 1/6, 1/3, 1/2, 2/3, etc.
  
  std::set<CGAL::Gmpq>  set;

  auto max_denom = 10;

  set.insert(1);

  for (auto denominator = 2; denominator <= max_denom ; ++denominator) 
  for (auto numerator = 1; numerator < denominator; ++numerator)
    set.insert(CGAL::Gmpq( numerator, denominator));

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

CGAL::Gmpq roundNice(const CGAL::Gmpq& x, const std::set<CGAL::Gmpq>& fracset, double eps) {

  // round the input to a "nice" fraction if it is eps-close (in absolute sense)
  // to such a fraction

  auto n = std::min(x.numerator().approximate_decimal_length(),
                    x.denominator().approximate_decimal_length());
  if (n < 20)
    return x;

  CGAL::Gmpz integral_part(x.to_double());

  auto frac = x - integral_part;

  auto res = find_nearest_nice_fraction(frac,fracset);

  if (res.second.to_double() < eps) {
    return integral_part + res.first;
  }
  return x;
}



CGAL::Gmpq roundQ(CGAL::Gmpq x, int precision) {

  // try to round the input to a fraction with precision as denominator
  //
  // (currently not used)
  auto n = std::min(x.numerator().approximate_decimal_length(),
                    x.denominator().approximate_decimal_length());
  if (n < 20)
    return x;

  auto z = (precision * x).to_double();
  auto rz = std::round(z);

  if (z - rz < 1e-11)
    return CGAL::Gmpq(rz, precision);
  return x;
}

template<class Kernel, class Container>
void print_polygon (const CGAL::Polygon_2<Kernel, Container>& P)
{
  typename CGAL::Polygon_2<Kernel, Container>::Vertex_const_iterator vit;
  std::cout << " " << P.size() << " vertices: {";
  for (vit = P.vertices_begin(); vit != P.vertices_end(); ++vit){
    if(vit != P.vertices_begin()) std::cout << ",";
    std::cout << "{" << (*vit).x() << ',' << (*vit).y() << "}";
  }
  std::cout << "} " << std::endl;
}

Polygon_2 convert_face_to_polygon(Face_handle f, int bound) {
  // Convert a voronoi region to a polygon of type Polygon_2
  // (We need this  because the CGAL::intersection function can, to the best of
  // our knowledge, not handle unbounded regions)
  //
  // If the region is unbounded, we clip it by adding an extra vertex on each
  // ray. The 'bound' parameter controls 'how far' the vertex will be placed

  Polygon_2 p;

  Ccb_halfedge_circulator ec_start = f->ccb();
  Ccb_halfedge_circulator ec = ec_start;

  do {
    if (ec->is_ray()) {
      Vector_2 v(ec->opposite()->face()->dual()->point(),
                 ec->face()->dual()->point());
      auto vp = v.perpendicular(CGAL::CLOCKWISE);

      auto t = bound;

      if (vp.squared_length() < 1){
        t *= int(1.0 / std::sqrt( CGAL::to_double(vp.squared_length())));
        // make sure that the extension will be not too short
      }

      if (ec->has_target()) {
        p.push_back(ec->target()->point() - t * vp);
        p.push_back(ec->target()->point());

      } else
        p.push_back(ec->source()->point() + t * vp);
    }
    else
      p.push_back(ec->target()->point());

  } while (++ec != ec_start);

  return p;
}

template<class Kernel, class Container>
void print_polygon_with_holes(const CGAL::Polygon_with_holes_2<Kernel, Container> & pwh)
{
  if (! pwh.is_unbounded()) {
    std::cout << "{ Outer boundary = "; 
    print_polygon (pwh.outer_boundary());
  } else
    std::cout << "{ Unbounded polygon." << std::endl;
  typename CGAL::Polygon_with_holes_2<Kernel,Container>::Hole_const_iterator hit;
  unsigned int k = 1;
  std::cout << " " << pwh.number_of_holes() << " holes:" << std::endl;
  for (hit = pwh.holes_begin(); hit != pwh.holes_end(); ++hit, ++k) {
    std::cout << " Hole #" << k << " = ";
    print_polygon (*hit);
  }
  std::cout << " }" << std::endl;
}

template <class Kernel, class Container>
Polygon_2
extract_poly(const CGAL::Polygon_with_holes_2<Kernel, Container> &poly) {
  // Extract a polygon from a more general object type
  // (CGAL::Polygon_with_holes_2), which is returned by several functions, like
  // minkowski sum. However, we never encounter polygons with holes.
  
  if (poly.is_unbounded())
    throw std::runtime_error("polygon is unbounded");

  if (poly.has_holes())
    throw std::runtime_error("polygon has holes");

  return poly.outer_boundary();
}

template <class Kernel, class Container>
Polygon_2 simple_intersect(const CGAL::Polygon_2<Kernel, Container> &PA,
                           const CGAL::Polygon_2<Kernel, Container> &PB) {
  // intersection of two overlapping convex polygons

  Pwh_list_2 intR;
  Pwh_list_2::const_iterator it;
  CGAL::intersection(PA, PB, std::back_inserter(intR));

  if (intR.size() != 1)
    throw std::runtime_error("intersection is disconnected");

  const auto &shape = *(intR.begin());

  return extract_poly(shape);
}

template <class Kernel, class Container>
Polygon_2 simple_union(const std::vector<CGAL::Polygon_2<Kernel, Container>> &input) {
  // union of a list of overlapping polygons

  std::list<Polygon_with_holes_2> _union_result;
  CGAL::join(input.cbegin(), input.cend(), std::back_inserter(_union_result));

  if (_union_result.size() != 1)
    throw std::runtime_error("union is disconnected");

  const auto &shape = *(_union_result.begin());
  return extract_poly(shape);

}

template <class Kernel, class Container>
Polygon_2 g(const CGAL::Polygon_2<Kernel, Container> &A,
       const CGAL::Polygon_2<Kernel, Container> &chS, const VD &voronoi) {
  // the main iteration

  // ignore the minkowski sum if the polygon is a singleton
  // (or if it has two vertices)
  Polygon_2 chS_plus_A;
  if(A.size() < 3)
   chS_plus_A = chS;
  else   
   chS_plus_A = extract_poly(CGAL::minkowski_sum_2(chS, A));


  // ensure that the union of the clipped voronoi regions cover chS+A

  auto max_abs_coord = [](const auto &point) {
    return CGAL::max(CGAL::abs(point.x()), CGAL::abs(point.y()));
  };
  auto it =
      std::max_element(chS_plus_A.vertices_begin(), chS_plus_A.vertices_end(),
                       [&max_abs_coord](const auto &a, const auto &b) {
        return (max_abs_coord(a) < max_abs_coord(b));
      });

  int bound = 100.0 * CGAL::to_double(max_abs_coord(*it));

  std::vector<Polygon_2> plist;

  // for each c, compute W_c = (chS + A) ∩ V(c) - c
  // where V(c) is the voronoi face belonging to delaunay-vertex ('site') c
  //
  // (however, we iterate over the faces and get c via the .dual() method)
  for (auto face_it = voronoi.faces_begin(); face_it != voronoi.faces_end();
       ++face_it) {

    auto c = face_it->dual()->point();
    Transformation translate_by_c(CGAL::TRANSLATION, Vector_2(c, Point_2(0, 0)));
    auto intersected_region = simple_intersect(chS_plus_A, convert_face_to_polygon(*face_it, bound));

    plist.push_back(transform(translate_by_c, intersected_region));
  }

  // take the union over the W_c sets for all c
  std::list<Polygon_with_holes_2> _union_result;
  CGAL::join(plist.cbegin(), plist.cend(), std::back_inserter(_union_result));

  if (_union_result.size() != 1)
    throw std::runtime_error("union is disconnected");

  const auto &shape = *(_union_result.begin());
  auto union_result = extract_poly(shape);

  // convex hull
  Polygon_2 convex_result;
  CGAL::ch_melkman( union_result.vertices_begin(), union_result.vertices_end(), std::back_inserter(convex_result));

  return convex_result;
}

template <class Kernel, class Container>
Polygon_2 round_vertices(const CGAL::Polygon_2<Kernel, Container> &P,int precision) {
  Polygon_2 res;
  for (auto i = P.vertices_begin(); i != P.vertices_end(); ++i) {
    auto p = *i;
    res.push_back(Point_2(roundQ(p.x(),precision), roundQ(p.y(),precision)));
  }
  return res;
}

template <class Kernel, class Container>
Polygon_2 round_vertices2(const CGAL::Polygon_2<Kernel, Container> &P,const std::set<CGAL::Gmpq>& fracset) {
  Polygon_2 res;
  for (auto i = P.vertices_begin(); i != P.vertices_end(); ++i) {
    auto p = *i;
    auto error = 1e-8;
    res.push_back(Point_2(roundNice(p.x(),fracset,error), roundNice(p.y(),fracset,error)));
  }
  return res;
}

template <class Kernel, class Container>
Polygon_2
remove_redundant_vertices(const CGAL::Polygon_2<Kernel, Container> &P) {
  // currently not used
  Polygon_2 reduced;
  auto i = P.vertices_circulator();
  auto start = i;
  do {
    if (!CGAL::collinear(*(i - 1), *i, *(i + 1)))
      reduced.push_back(*i);
  } while (++i != start);
  return reduced;
}

CGAL::Gmpq get_frac(const nlohmann::json &j) {
  if (j.is_string())
    return j.get<std::string>();
  else
    return j.get<int>();
}

auto precompute_convex_hulls_and_voronoi_diagrams(const nlohmann::json &user_input)
{

  std::vector<std::pair<Polygon_2,VD>> cvs; // convex hulls and voronoi diagrams 

  for (auto& point_set : user_input["point_sets"])
  {

    VD vd;
    for (const auto &coord_pair : point_set) {
      Point_2 p(get_frac(coord_pair[0]), get_frac(coord_pair[1]));
      vd.insert(p);
    }

    assert(vd.is_valid());
    Polygon_2 chS;
    CGAL::ch_graham_andrew(vd.sites_begin(), vd.sites_end(),
                           std::back_inserter(chS));

    cvs.emplace_back(std::make_pair(chS,vd));
  }

  return cvs;
}


void echo_parsed_input(const std::vector<std::pair<Polygon_2,VD>>& cvs)
{
  auto sz = cvs.size();
  std::cout << "Read " << sz << " point " << ((sz==1) ? "set" : "sets") << ":" << std::endl;

  for (const auto& pair : cvs)
  {
    auto vd = pair.second;
    Polygon_2 P;
    for (auto it = vd.sites_begin(); it != vd.sites_end(); ++it)
      P.push_back(*it);
    print_polygon(P);
  }
}

Polygon_2 parse_or_create_initial_set(const nlohmann::json &user_input)
{
  Polygon_2 A;
  if (user_input.count("A") == 1){
    for (const auto &coord_pair :user_input["A"])
    {
      Point_2 p(get_frac(coord_pair[0]), get_frac(coord_pair[1]));
      A.push_back(p);
    }
  }
  else
    A.push_back(Point_2(0, 0));

  return A;
}

int main(int argc, char* argv[])
{
  using json = nlohmann::json;

  if(argc < 2)
  {
    std::cerr << "Usage: " << argv[0] << " <points file>" << std::endl;
    return 1;
  }

  // read points and build Voronoi diagram
  std::ifstream ifs(argv[1]);
  assert( ifs );

  json j;
  ifs >> j;
  ifs.close();

  auto cvs = precompute_convex_hulls_and_voronoi_diagrams(j); // convex hulls and voronoi diagrams 
  echo_parsed_input(cvs);
  
  auto A = parse_or_create_initial_set(j);
  std::cout << "We will use the following polygon as starting set." << std::endl;
  print_polygon(A);

  Polygon_2 Aprev;

  // needed for rounding
  auto s = build_nice_fraction_list();

  // run fixed point
  std::cout << "Starting fixed point iteration..." << std::endl;
  auto iter = 0; 
  while (A != Aprev) {
    Aprev = A;

    std::vector<Polygon_2> plist;
    for (auto &p : cvs) {

      plist.push_back(g(A, p.first, p.second));
    }

    Polygon_2 union_result;
    union_result = round_vertices2(simple_union(plist), s);

    Polygon_2 convex_result;
    CGAL::ch_melkman(union_result.vertices_begin(), union_result.vertices_end(),
                     std::back_inserter(convex_result));

    A = convex_result;
    iter++;
    print_polygon(A);
  }

  std::cout << "Found invariant set after " << iter
            << ((iter == 1) ? " iteration." : " iterations.") << std::endl;

  return 0;
}


