/*
 * Copyright Ahnaf Siddiqui and Sameer Varma.
 */

typedef dtreal double; // test performance and accuracy of double vs float


struct dTriangulation {
	dtreal *points; // coordinates of input points (2 ordered reals, x and y, per point)
	int npoints;

	int *triangles; // list of delaunay triangles as groups of three point indexes 
					// (index = order of point in given points array)
	int ntriangles;
};


void dtriangulate(struct dTriangulation *tri);
/* Triangulates the points given in tri using Delaunay triangulation,
 * storing the resulting triangles in tri (see struct above).
 * Memory is allocated for tri->triangles.
 */
