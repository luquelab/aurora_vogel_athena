# What
This folder contains source code files that require compilation.

Solver.py contains my attempt at using sympy solvers to solve the
system of 27 quadratics with 27 variables.
-> Comments are above each attempt 
   to indicate results
-> I worked through this stuff via the interactive
   python window in MS Visual Studio. 
   I attempted to wrap things up into some functions for easier use
   (depends on how you like to work).
-> The triangle distances were hardcoded (if the easiest case can't be solved,
   no sense in generalizing it).

minimizer.py contains my attempt at minimizing the distance function which
has the distances b/w the actual verts minus there theoretical
distances
-> I tried a couple minimization methods in the minimize function
   but results were the same
-> feeding in even a slightly perturbed set of verts
   (only one vertex slightly off) was not solvable
   (perhaps it's at a local minimum...?).

folding.py can be ignored (I thought I could model the folding process
however it's more complicated than rotating about two axes at a time).