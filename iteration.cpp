// standard includes
#include <iostream>
#include <cassert>
#include <cmath>
#include <cctype>
#include <stdexcept>
#include <algorithm>
#include <functional>
#include <unordered_map>
#include <string>

#include <CGAL/Aff_transformation_2.h>
#include <CGAL/ch_graham_andrew.h>
#include <CGAL/ch_melkman.h>

#include "json.hpp" // nlohmann's json lib
#include "polygon.hpp"
#include "rounding.hpp"

Polygon_2 g_inner(const Polygon_2 &chS_plus_A, const VD &voronoi)
{
  // ensure that the union of the clipped voronoi regions cover chS+A
  int bound = 100.0 * find_bounding_radius(chS_plus_A);

  // for each c, compute W_c = (chS + A) ∩ V(c) - c
  // where V(c) is the voronoi face belonging to delaunay-vertex ('site') c
  //
  // (however, we iterate over the faces and get c via the .dual() method)
  std::vector<Polygon_2> plist;
  for (auto face_it = voronoi.faces_begin(); face_it != voronoi.faces_end();
       ++face_it) {
    auto c = face_it->dual()->point();
    CGAL::Aff_transformation_2<K> translate_by_c(CGAL::TRANSLATION, CGAL::Vector_2<K>(c, Point_2(0, 0)));
    auto intersected_region = intersection(chS_plus_A, convert_face_to_polygon(*face_it, bound));
    plist.push_back(transform(translate_by_c, intersected_region));
  }

  // take the union over the W_c sets for all c
  return union_of_overlapping_polygons(plist);
}

Polygon_2 g(const Polygon_2 &A, const Polygon_2 &chS, const VD &voronoi) {
  //non convex G-iteration
  Polygon_2 chS_plus_A = minkowski_sum_between_nonconv_and_conv(A, chS);
  return g_inner(chS_plus_A, voronoi);
}

Polygon_2 G(const Polygon_2 &A, const Polygon_2 &chS, const VD &voronoi) {
  // convex G-iteration
  //
  // assumes A is convex
  Polygon_2 U = g_inner(minkowski_sum_between_convex_sets(chS, A), voronoi);
  Polygon_2 convex_result;
  CGAL::ch_melkman(U.vertices_begin(), U.vertices_end(), std::back_inserter(convex_result));
  return convex_result;
}

Polygon_2 f(const Polygon_2 &D, const Polygon_2 &chS, const VD &voronoi) {
  // non-convex f-iteration
  return minkowski_sum_between_nonconv_and_conv(g_inner(D, voronoi), chS);
}

Polygon_2 F(const Polygon_2 &D, const Polygon_2 &chS, const VD &voronoi) {
  // convex F-iteration
  Polygon_2 chS_plus_UD = minkowski_sum_between_nonconv_and_conv(g_inner(D, voronoi), chS);
  Polygon_2 convex_result;
  CGAL::ch_melkman(chS_plus_UD.vertices_begin(), chS_plus_UD.vertices_end(), std::back_inserter(convex_result));
  return convex_result;
}

CGAL::Gmpq get_frac(const nlohmann::json &j) {
  if (j.is_string())
    return j.get<std::string>();
  else
    return j.get<int>();
}

auto precompute_convex_hulls_and_voronoi_diagrams(
    const nlohmann::json &user_input) {

  std::vector<std::pair<Polygon_2, VD>> cvs;
  // convex hulls and voronoi diagrams

  for (auto &point_set : user_input["point_sets"]) {

    VD vd;
    for (const auto &coord_pair : point_set) {
      Point_2 p(get_frac(coord_pair[0]), get_frac(coord_pair[1]));
      vd.insert(p);
    }

    assert(vd.is_valid());
    Polygon_2 chS;
    CGAL::ch_graham_andrew(vd.sites_begin(), vd.sites_end(),
                           std::back_inserter(chS));

    cvs.emplace_back(std::make_pair(chS, vd));
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

Polygon_2
parse_or_create_initial_set(const nlohmann::json &user_input,
                            const std::string &iter,
                            const std::vector<std::pair<Polygon_2, VD>> &cvs) {
  Polygon_2 A;
  if (user_input.count("A") == 1){
    for (const auto &coord_pair :user_input["A"])
    {
      Point_2 p(get_frac(coord_pair[0]), get_frac(coord_pair[1]));
      A.push_back(p);
    }
  }
  else if (std::toupper(iter[0]) == 'G')
    A.push_back(Point_2(0, 0));

  else if (std::toupper(iter[0]) == 'F')
    A = cvs[0].first;
  
  else throw std::runtime_error("ERROR: failed to create initial set");

  return A;
}

template <template <typename...> class Iterable>
void run_fixed_point_iteration(
    std::function<Polygon_2(const Polygon_2 &, const Polygon_2 &, const VD &)>& fun,
    Polygon_2 &A, const Iterable<std::pair<Polygon_2, VD>> &cvs,
    int digit_threshold, double error, bool apply_conv_hull) {

  // needed for rounding
  auto s = build_nice_fraction_list(1000);

  std::vector<int> rounding_bookkeep;
  Polygon_2 A_original;
  std::cout << "Starting fixed point iteration..." << std::endl;
  auto iter = 0; 
  while (true) {
    A_original = A;

    std::vector<Polygon_2> plist;
    for (auto &p : cvs) {
      plist.push_back(fun(A, p.first, p.second));
    }
    A = union_of_overlapping_polygons(plist);

    if (A == A_original)
      // did we converge?
      break;

    if (iter % 5 == 0) { // do not round in every iteration
      auto V = round_vertices(A, s, digit_threshold, error);

      A = V.first;

      if (V.second) {
        std::cout << "Applied rounding in iteration " << iter << std::endl;
        A = remove_redundant_vertices(A);
        rounding_bookkeep.push_back(iter);
      }
    }

    iter++;
    if (iter % 1 == 0) {
      std::cout << "iteration: " << iter << std::endl;
      print_polygon(A);
    }

  }

  std::cout << "Found invariant set after " << iter
            << ((iter == 1) ? " iteration." : " iterations.") << std::endl;

  if (rounding_bookkeep.size() == 0)
    std::cout << "Coefficient rounding was not necessary and has not been applied.\n";
  else
  {
    std::cout << "Coefficient rounding has been applied in the following " << rounding_bookkeep.size() << " iterations:\n[";
    for (int i : rounding_bookkeep)
      std::cout << i << ", ";
    std::cout << "]\n\n";
  }

  std::cout << "The coordinates of the invariant polygon are:" << std::endl;
  print_polygon(A);
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

  double error(1e-8);
  if (j.count("rounding_error") == 1) {
    error = j["rounding_error"].get<double>();
  }

  int digit_threshold(10);
  if (j.count("digit_threshold") == 1) {
    digit_threshold = j["digit_threshold"].get<int>();
  }

  std::unordered_map<std::string,
                     std::function<Polygon_2(const Polygon_2 &,
                                             const Polygon_2 &, const VD &)>>
      function_map = {{"g", g}, {"G", G}, {"f", f}, {"F", F}};

  std::string iter_name("g");
  if (j.count("iteration") == 1) {
    iter_name = j["iteration"].get<std::string>();

    // validate input
    if (function_map.find(iter_name) == function_map.end())
    {
      std::cout << "ERROR: Unknown user-specified iteration: " << iter_name << "\n";
      std::cout << "Available iteration types:\n"; 
      for(const auto& kv: function_map)
        std::cout << "    " << kv.first << '\n';
      return -1;
    }
  }

  auto A = parse_or_create_initial_set(j, iter_name, cvs);
  std::cout << "We will use the following polygon as starting set." << std::endl;
  print_polygon(A);

  std::cout << "Performing the " << iter_name << "-iteration.\n";

  bool apply_conv_hull = (iter_name[0] == 'F');
  run_fixed_point_iteration(function_map[iter_name], A, cvs, digit_threshold, error, apply_conv_hull);

  return 0;
}


