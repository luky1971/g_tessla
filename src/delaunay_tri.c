/*
 * Copyright Ahnaf Siddiqui and Sameer Varma.
 *
 * Delaunay triangulation implementation based on the algorithms presented in
 * "Two Algorithms for Constructing the Delaunay Triangulation," 
 	International Journal of Computer and Information Science 9(3):219-242, 1980
 * and
 * "Primitives for the Manipulation of General Subdivisions and the Computation of Voronoi Diagrams," 
 	ACM Transactions on Graphics 4(2):74-123, April 1985
 *
 */

#include "delaunay_tri.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define DTEPSILON 1e-8 // explore different epsilon values


struct vert {
	int index;
	struct vertNode *adj;
};

struct vertNode {
	struct vert *v;
	struct vertNode *prev, *next;
};


int comparePoints(const void *a, const void *b);
static void ord_dtriangulate(struct vert *v, int ia, int ib);


static struct dTriangulation *mtri;


// lexicographically compares the points in mtri->points given by indexes a and b.
int comparePoints(const void *a, const void *b) {
	int ia = *((int*)a), ib = *((int*)b);
	dtreal xa = mtri->points[ia*2];
	dtreal ya = mtri->points[ia*2+1];
	dtreal xb = mtri->points[ib*2];
	dtreal yb = mtri->points[ib*2+1];

	dtreal diff = xa - xb;

	if(diff < DTEPSILON && diff > -DTEPSILON) {
		diff = ya - yb;
	}

	return (1 - 2 * signbit(diff));
}

void dtriangulate(struct dTriangulation *tri) {
	mtri = tri;
	int *ind = (int*)malloc(mtri->npoints * sizeof(int));

	for(int i = 0; i < mtri->npoints; ++i) {
		ind[i] = i;
	}

	// sort point indexes lexicographically by point coordinates
	qsort(ind, mtri->npoints, sizeof(int), comparePoints);

	// construct vertex structures of points in sorted order
	struct vert *v = (struct vert*)malloc(mtri->npoints * sizeof(struct vert));

	for(int i = 0; i < mtri->npoints; ++i) {
		v[i].index = ind[i];
		v[i].adj = NULL;
	}

	free(ind);

	// DEBUG
	// FILE *f = fopen("points.txt", "w");

	// for(int i = 0; i < mtri->npoints; ++i) {
	// 	fprintf(f, "%f, %f\n", mtri->points[2 * v[i].index], mtri->points[2 * v[i].index + 1]);
	// }

	// fclose(f);
	//

	// triangulate the sorted points
	ord_dtriangulate(v, 0, mtri->npoints - 1);

	// TODO: free vertnodes in adjacency lists!
	free(v);
}

static void ord_dtriangulate(struct vert *v, int ia, int ib) {
	if(ia >= ib)
		return;
	if(ib - ia == 1) {
		// num points = 2. Handle this base case
	}
	else if(ib - ia == 2) {
		// num points = 3. Handle this base case
	}
	else {
		int mid = (ia + ib) / 2;
		ord_dtriangulate(v, ia, mid);
		ord_dtriangulate(v, mid + 1, ib);
	}
}
