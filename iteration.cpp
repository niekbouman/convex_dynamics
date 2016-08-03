// standard includes
#include <iostream>
#include <fstream>
#include <cassert>
#include <cmath>
#include <stdexcept>
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

//#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/ch_graham_andrew.h>
#include <CGAL/ch_melkman.h>
#include <CGAL/Aff_transformation_2.h>

// typedefs for defining the adaptor
// typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Exact_predicates_exact_constructions_kernel K;
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
typedef CGAL::Polygon_with_holes_2<K> Polygon_with_holes_2;
typedef std::list<Polygon_with_holes_2> Pwh_list_2;

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

Polygon_2 convert_face_to_polygon(Face_handle *f, int bound) {
  // Convert a voronoi region to a polygon of type Polygon_2
  // (We need this  because the CGAL::intersection function can, to the best of
  // our knowledge, not handle unbounded regions)
  //
  // If the region is unbounded, we clip it by adding an extra vertex on each
  // ray. The 'bound' parameter controls 'how far' the vertex will be placed

  Polygon_2 p;

  Ccb_halfedge_circulator ec_start = (*f)->ccb();
  Ccb_halfedge_circulator ec = ec_start;

  do {
    if (ec->is_ray()) {
      Vector_2 v(ec->opposite()->face()->dual()->point(),
                 ec->face()->dual()->point());
      auto vp = v.perpendicular(CGAL::CLOCKWISE);

      if (vp.squared_length() < 1){
        bound *= int(1.0 / std::sqrt( CGAL::to_double(vp.squared_length())));
        // make sure that the extension will be  not too short
      }

      if (ec->has_target()) {
        p.push_back(ec->target()->point() - bound * vp);
        p.push_back(ec->target()->point());

      } else
        p.push_back(ec->source()->point() + bound * vp);
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

  std::vector<Polygon_2> plist;

  // for each c, compute W_c = (chS + A) ∩ V(c) - c
  for (auto c = voronoi.sites_begin(); c != voronoi.sites_end(); ++c ) {

    Transformation translate_by_c(CGAL::TRANSLATION, Vector_2(*c,Point_2(0,0)));
    auto result = voronoi.locate(*c);
    if (Face_handle *f = boost::get<Face_handle>(&result)) {

      auto intersected_region =
          simple_intersect(chS_plus_A, convert_face_to_polygon(f, 1000));
      plist.push_back(transform(translate_by_c, intersected_region));
    }
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


int main()
{

  // read points and build Voronoi diagram
  std::ifstream ifs("points.cin");
  assert( ifs );
  VD vd;
  Site_2 t;
  while ( ifs >> t ) { vd.insert(t); }
  ifs.close();
  assert( vd.is_valid() );

  // compute convex hull of the points
  Polygon_2 chS;
  CGAL::ch_graham_andrew(vd.sites_begin(),vd.sites_end(), std::back_inserter(chS));

  Polygon_2 A;
  A.push_back(Point_2(0,0));
  Polygon_2 Aprev;

  // run fixed point
  while(A != Aprev) {
    Aprev = A;
    A = g(A,chS,vd);
    print_polygon(A);
  } 

  return 0;
}
