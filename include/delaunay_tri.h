/*
 * Copyright Ahnaf Siddiqui and Sameer Varma.
 *
 * Delaunay triangulation implementation based on the algorithms presented by
 *
 * 	Der-Tsai Lee and Bruce J. Schachter, 
 	"Two Algorithms for Constructing the Delaunay Triangulation," 
	International Journal of Computer and Information Science 9(3):219-242, 1980
 *
 * and
 *
 * 	Leonidas J. Guibas and Jorge Stolfi, 
	"Primitives for the Manipulation of General Subdivisions and the Computation of Voronoi Diagrams," 
	ACM Transactions on Graphics 4(2):74-123, April 1985
 *
 */

typedef double dtreal; // test performance and accuracy of double vs float


struct dTriangulation {
	dtreal *points; // coordinates of input points (2 ordered reals, x and y, per point)
	int npoints;

	int *triangles; // list of delaunay triangles as groups of three point indexes 
					// (index = order of point in given points array)
	int ntriangles;
	int nverts; // equivalent to the number of non-duplicate input points
};


void dtriangulate(struct dTriangulation *tri);
/* Triangulates the points given in tri using Delaunay triangulation,
 * storing the resulting triangles in tri (see struct above).
 * Memory is allocated for tri->triangles.
 */
