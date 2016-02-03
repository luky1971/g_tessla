/*
 * Copyright Ahnaf Siddiqui and Sameer Varma.
 */

#include "delaunay_tri.h"

#include <math.h>
#include <stdlib.h>

#define DTEPSILON 1e-8 // explore different epsilon values

int comparePoints(const void *a, const void *b);
static void ord_dtriangulate();


static struct dTriangulation *mtri;


// lexicographically compares the points in mtri->points given by indexes a and b.
int comparePoints(const void *a, const void *b) {
	dtreal xa = mtri->points[*((int*)a) * 2];
	dtreal ya = mtri->points[*((int*)a) * 2 + 1];
	dtreal xb = mtri->points[*((int*)b) * 2];
	dtreal yb = mtri->points[*((int*)a) * 2 + 1];

	dtreal diff = xa - xb;

	if(diff < DTEPSILON && diff > -DTEPSILON) {
		diff = ya - yb;
	}

	return -1 * signbit(diff);
}

void dtriangulate(struct dTriangulation *tri) {
	mtri = tri;
	// int *ind = (int*)malloc(mtri->npoints * sizeof(int));

	// for(int i = 0; i < mtri->npoints; ++i) {
	// 	ind[i] = i;
	// }

	// sort point indexes lexicographically by point coordinates
	// qsort(ind, mtri->npoints, sizeof(int), comparePoints);

	// free(ind);

	
}
