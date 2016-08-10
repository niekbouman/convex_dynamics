# convex_dynamics
Computation of the minimal convex invariant error set for the error-diffusion algorithm in 2D, given a collection of sets of points.

Uses the CGAL C++ library.

Point set can be given in a JSON file with the collection of point sets as a triply-nested array. 
Example: below we have a collection of three point sets:

   [
   [[-1,1],[1,1],[1,-1],[-1,-1],[-2,-1]],
   [[-1,1],[1,1],[1,-1],[-1,-1],[-2,-1],[-2,5]],
   [[-2,2],[1,5],[1,-1]]
   ]

Save this in a file, and run `./iter <filename>`
