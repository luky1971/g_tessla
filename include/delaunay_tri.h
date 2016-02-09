/*
 * Copyright 2016 Ahnaf Siddiqui and Sameer Varma
 *
 * Delaunay triangulation implementation based on the algorithms presented by
 *
 * Lee, D.T. and Schachter, B.J. 
 * Two algorithms for constructing a Delaunay triangulation. 
 * International Journal of Computer & Information Sciences 1980;9(3):219-242.
 *
 * and
 *
 * Guibas, L. and Stolfi, J. 
 * Primitives for the manipulation of general subdivisions and the computation of Voronoi. 
 * ACM Trans. Graph. 1985;4(2):74-123.
 *
 * This implementation uses exact arithmetic routines and geometric predicates provided by
 *
 * Shewchuk, J.R. 1996. 
 * Routines for Arbitrary Precision Floating-point Arithmetic and Fast Robust Geometric Predicates.
 */

#include "predicates.h"

typedef REAL dtreal; // test performance and accuracy of double vs float


struct dTriangulation {
	dtreal *points; // coordinates of input points (2 ordered reals, x and y, per point)
	int npoints;

	int *triangles; // list of delaunay triangles as groups of three point indexes 
					// (index = order of point in given points array)
	int ntriangles;
	int nverts; // equivalent to the number of non-duplicate input points
};


void dtinit();
/* Call this once before calling dtriangulate()
 */

void dtriangulate(struct dTriangulation *tri);
/* Triangulates the points given in tri using Delaunay triangulation,
 * storing the resulting triangles in tri (see struct above).
 * Memory is allocated for tri->triangles.
 */
