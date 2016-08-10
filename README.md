# convex_dynamics
Computation of the minimal convex invariant error set for the error-diffusion algorithm in 2D, given a collection of sets of points.

Uses the CGAL C++ library.

Point set can be given in a JSON object with the collection of point sets with key "point_sets" as an array of arrays of pairs of coordinates. Coordinates can be integers, or fractions, denoted as strings, like "4/5".
Example: below we have a collection of three point sets:

```
{"point_sets":[
  [[-1,1],[1,1],[1,"-3/8"],[-1,-1],[-2,-1]],
  [[-1,1],[1,"1/2],[1,-1],[-1,-1],[-2,-1],[-2,5]],
  [[-2,2],[1,5],[1,-1]]
]}
```
Save this in a file, and run `./iter <filename>`

Optionally, the initial set can be specified with key "A":
```
{
 "A":[[1,1],[-1,1],[-1,-1],[1,-1]],
 "point_sets":[...] 
}
```
If "A" is omitted, the singleton set [0,0] is used.
