#include "polygon.hpp"

//#include <CGAL/Direction_2.h>

#include <CGAL/Polygon_nop_decomposition_2.h>
#include <CGAL/Small_side_angle_bisector_decomposition_2.h>
#include <CGAL/Minkowski_sum_2.h>

// typedef for the result type of the point location

//using Site_2 = AT::Site_2;
//using Locate_result = VD::Locate_result;
//using Vertex_handle = VD::Vertex_handle;
//using Halfedge_handle = VD::Halfedge_handle;



void print_polygon(const Polygon_2 &P) {
  using std::cout;
  cout << " " << P.size() << " vertices: {";
  for (auto vit = P.vertices_begin(); vit != P.vertices_end(); ++vit) {
    if (vit != P.vertices_begin())
      cout << ",";
    cout << "{" << (*vit).x() << ',' << (*vit).y() << "}";
  }
  cout << "}\n" ;
}

Polygon_2 convert_face_to_polygon(Face_handle f, int bound) {
  // Convert a voronoi region to a polygon of type Polygon_2
  // (We need this  because the CGAL::intersection function can, to the best of
  // our knowledge, not handle unbounded regions)
  //
  // If the region is unbounded, we clip it by adding an extra vertex on each
  // ray. The 'bound' parameter controls 'how far' the vertex will be placed

  Polygon_2 p;

  using Ccb_halfedge_circulator = VD::Ccb_halfedge_circulator;
  Ccb_halfedge_circulator ec_start = f->ccb();
  Ccb_halfedge_circulator ec = ec_start;

  do {
    if (ec->is_ray()) {

      CGAL::Vector_2<K> v(ec->opposite()->face()->dual()->point(),
                          ec->face()->dual()->point());
      auto vp = v.perpendicular(CGAL::CLOCKWISE);

      auto t = bound;

      if (vp.squared_length() < 1) {
        t *= int(1.0 / std::sqrt(CGAL::to_double(vp.squared_length())));
        // make sure that the extension will be not too short
      }

      if (ec->has_target()) {
        p.push_back(ec->target()->point() - t * vp);
        p.push_back(ec->target()->point());

      } else
        p.push_back(ec->source()->point() + t * vp);
    } else
      p.push_back(ec->target()->point());

  } while (++ec != ec_start);

  return p;
}

/*
void print_polygon_with_holes(const Polygon_with_holes_2& pwh)
{
  if (! pwh.is_unbounded()) {
    std::cout << "{ Outer boundary = "; 
    print_polygon (pwh.outer_boundary());
  } else
    std::cout << "{ Unbounded polygon." << std::endl;
  unsigned int k = 1;
  std::cout << " " << pwh.number_of_holes() << " holes:" << std::endl;
  for (auto hit = pwh.holes_begin(); hit != pwh.holes_end(); ++hit, ++k) {
    std::cout << " Hole #" << k << " = ";
    print_polygon (*hit);
  }
  std::cout << " }" << std::endl;
}
*/

Polygon_2 extract_poly(const Polygon_with_holes_2 &poly) {
  // Extract a polygon from a more general object type
  // (CGAL::Polygon_with_holes_2), which is returned by several functions, like
  // minkowski sum. However, we never encounter polygons with holes.
  
  if (poly.is_unbounded())
    throw std::runtime_error("polygon is unbounded");

  if (poly.has_holes())
    throw std::runtime_error("polygon has holes");

  return poly.outer_boundary();
}

Polygon_2 intersection(const Polygon_2 &PA, const Polygon_2 &PB) {
  // intersection of two overlapping polygons

  std::list<Polygon_with_holes_2> intR;
  CGAL::intersection(PA, PB, std::back_inserter(intR));

  if (intR.size() != 1)
    throw std::runtime_error("intersection is disconnected");

  const auto &shape = *(intR.begin());

  return extract_poly(shape);
}

Polygon_2 union_of_overlapping_polygons(const std::vector<Polygon_2> &input) {
  // union of a list of overlapping polygons

  std::list<Polygon_with_holes_2> _union_result;
  CGAL::join(input.cbegin(), input.cend(), std::back_inserter(_union_result));

  if (_union_result.size() != 1)
    throw std::runtime_error("union is disconnected");

  const auto &shape = *(_union_result.begin());
  return remove_collinear_vertices(extract_poly(shape));

}

Polygon_2 remove_collinear_vertices(const Polygon_2 &P) {
  Polygon_2 reduced;
  auto i = P.vertices_circulator();
  auto start = i;
  do {
    if (!CGAL::collinear(*(i - 1), *i, *(i + 1)))
      reduced.push_back(*i);
  } while (++i != start);
  return reduced;
}

Polygon_2 remove_repeated_vertices(const Polygon_2 &P) {
  Polygon_2 reduced;
  auto i = P.vertices_circulator();
  auto start = i;
  do {
    if (*i != *(i + 1))
      reduced.push_back(*i);
  } while (++i != start);
  return reduced;
}

Polygon_2 remove_redundant_vertices(Polygon_2 P) {
  Polygon_2 Pp;
  do {
    Pp = P;
    P = remove_collinear_vertices(P);
    P = remove_repeated_vertices(P);
  } while (Pp != P);
  return P;
}

Polygon_2 minkowski_sum_between_nonconv_and_conv(const Polygon_2 &P,
                                                 const Polygon_2 &Q)
{
  CGAL::Small_side_angle_bisector_decomposition_2<K> ssab_decomp;
  CGAL::Polygon_nop_decomposition_2<K> nop;
  return extract_poly(CGAL::minkowski_sum_2(P, Q, ssab_decomp, nop));
}

Polygon_2 minkowski_sum_between_convex_sets(const Polygon_2 &P,
                                            const Polygon_2 &Q) {
  CGAL::Polygon_nop_decomposition_2<K> nop;
  return extract_poly(CGAL::minkowski_sum_2(P, Q, nop));
}

double find_bounding_radius(const Polygon_2 &P) {

  auto max_abs_coord = [](const auto &point) {
    return CGAL::max(CGAL::abs(point.x()), CGAL::abs(point.y()));
  };
  auto it = std::max_element(P.vertices_begin(), P.vertices_end(),
                             [&max_abs_coord](const auto &a, const auto &b) {
                               return (max_abs_coord(a) < max_abs_coord(b));
                             });

  return CGAL::to_double(max_abs_coord(*it));
}

