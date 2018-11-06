#ifndef POLYGONHPP
#define POLYGONHPP

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Gmpq.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Voronoi_diagram_2.h>
#include <CGAL/Delaunay_triangulation_adaptation_traits_2.h>
#include <CGAL/Delaunay_triangulation_adaptation_policies_2.h>

using K = CGAL::Simple_cartesian<CGAL::Gmpq>;
using Polygon_2 = CGAL::Polygon_2<K>;
using Polygon_with_holes_2 = CGAL::Polygon_with_holes_2<K>;

using DT = CGAL::Delaunay_triangulation_2<K>;
using AT = CGAL::Delaunay_triangulation_adaptation_traits_2<DT>;
using AP = CGAL::Delaunay_triangulation_caching_degeneracy_removal_policy_2<DT>;
using VD = CGAL::Voronoi_diagram_2<DT, AT, AP>;

using Face_handle = VD::Face_handle;
using Point_2 = AT::Point_2;

void print_polygon(const Polygon_2 &P);

Polygon_2 convert_face_to_polygon(Face_handle f, int bound);
// Convert a voronoi region to a polygon of type Polygon_2
// (We need this  because the CGAL::intersection function can, to the best of
// our knowledge, not handle unbounded regions)
//
// If the region is unbounded, we clip it by adding an extra vertex on each
// ray. The 'bound' parameter controls 'how far' the vertex will be placed

Polygon_2 extract_poly(const Polygon_with_holes_2 &poly);
// Extract a polygon from a more general object type
// (CGAL::Polygon_with_holes_2), which is returned by several functions, like
// minkowski sum. However, we never encounter polygons with holes.

Polygon_2 intersection(const Polygon_2 &PA, const Polygon_2 &PB);
// intersection of two overlapping polygons

Polygon_2 union_of_overlapping_polygons(const std::vector<Polygon_2> &input);
// union of a list of overlapping polygons

Polygon_2 remove_collinear_vertices(const Polygon_2 &P);
Polygon_2 remove_repeated_vertices(const Polygon_2 &P);
Polygon_2 remove_redundant_vertices(Polygon_2 P);

Polygon_2 minkowski_sum_between_nonconv_and_conv(const Polygon_2 &P,
                                                 const Polygon_2 &Q);
Polygon_2 minkowski_sum_between_convex_sets(const Polygon_2 &P,
                                            const Polygon_2 &Q);

double find_bounding_radius(const Polygon_2 &P);


#endif
