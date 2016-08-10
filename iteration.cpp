// standard includes
#include <iostream>
#include <utility>

#include <fstream>
#include <cassert>
#include <cmath>
#include <stdexcept>
#include <algorithm>
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

CGAL::Gmpq roundQ(CGAL::Gmpq x, int precision) {

  // try to round the input to a fraction with precision as denominator
  
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
    std::cout << "{" << (*vit).x().to_double() << ',' << (*vit).y().to_double() << "}";
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
Polygon_2
remove_redundant_vertices(const CGAL::Polygon_2<Kernel, Container> &P) {
  Polygon_2 reduced;
  auto i = P.vertices_circulator();
  auto start = i;
  do {
    if (!CGAL::collinear(*(i - 1), *i, *(i + 1)))
      reduced.push_back(*i);
  } while (++i != start);
  return reduced;
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

  std::cout << "Read " << j.size() << " point set(s) from "<< argv[1] << "." << std::endl;

 
  std::vector<std::pair<Polygon_2,VD>> cvs; // convex hulls and voronoi diagrams 
  for (auto& point_set : j)
  {

    VD vd;
    for (const auto &coord_pair : point_set) {
      Point_2 p(coord_pair[0].get<int>(), coord_pair[1].get<int>());
      vd.insert(p);
    }

    assert(vd.is_valid());
    Polygon_2 chS;
    CGAL::ch_graham_andrew(vd.sites_begin(), vd.sites_end(),
                           std::back_inserter(chS));


    cvs.emplace_back(std::make_pair(chS,vd));
  }

  Polygon_2 A;
  A.push_back(Point_2(0, 0));
  Polygon_2 Aprev;

  // run fixed point
  while (A != Aprev) {
    Aprev = A;

    std::vector<Polygon_2> plist;
    for (auto &p : cvs) {

      plist.push_back(g(A, p.first, p.second));
    }

    Polygon_2 union_result;
    union_result = round_vertices(simple_union(plist), 100);

    Polygon_2 convex_result;
    CGAL::ch_melkman(union_result.vertices_begin(), union_result.vertices_end(),
                     std::back_inserter(convex_result));

    A = convex_result;
    print_polygon(A);
  }

  return 0;
}


