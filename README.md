# convex_dynamics
Computation of the minimal convex invariant error set for the error-diffusion algorithm in 2D, given a collection of sets of points.

Uses the CGAL C++ library.

Point set can be given in a JSON object with the collection of point sets with key "point_sets" as an array of arrays of pairs of coordinates. Coordinates can be integers, or fractions, denoted as strings, like "4/5".
Example: below we have a collection of three point sets:

```
{"point_sets":[
  [[-1,1],[1,1],[1,"-3/8"],[-1,-1],[-2,-1]],
  [[-1,1],[1,"1/2"],[1,-1],[-1,-1],[-2,-1],[-2,5]],
  [[-2,2],[1,5],[1,-1]]
]}
```
Save this in a file, and run `./iter <filename>`

## Options
The following options can be specified:
 - "iteration" : specify which iteration to use:
   - "g" : the G-iteration (default)
   - "G" : the convexified G-iteration
   - "f" : the F-iteration
   - "F" : the convexified F-iteration
 - "A" : the initial set to use 
   For example,
   ```
   {
     "A":[[1,1],[-1,1],[-1,-1],[1,-1]],
     "point_sets":[...] 
   }
   ```
   If "A" is omitted, the singleton set [0,0] is used for g and G iteration, 
   while the convex hull of the first point set is used for the f and F iterations.

 - "digit_threshold" : threshold specifying how many digits the numerator or denominator of a vertex coordinate should have in order to enable coefficient rounding 
 - "rounding_error" : round a coordinate when it is within rounding_error from a small fraction  

