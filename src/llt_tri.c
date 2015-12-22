/*
 * Copyright Ahnaf Siddiqui and Sameer Varma.
 */

#include "llt_tri.h"

#include <float.h>
#include <string.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "gkut_io.h"
#include "gkut_log.h"
#include "smalloc.h"


void llt_tri_area(const char *traj_fname, const char *ndx_fname, output_env_t *oenv, 
	real **areas, int *nframes, int *natoms, unsigned char flags) {
	rvec **pre_x, **x;
	real area;

	read_traj(traj_fname, &pre_x, nframes, natoms, oenv);

	// Filter trajectory by index file if present
	if(ndx_fname != NULL) {
		ndx_filter_traj(ndx_fname, &pre_x, &x, *nframes, natoms);
	}
	else {
		x = pre_x;
	}

	// Calculate triangulated surface area for every frame
	snew(*areas, *nframes);
	if(flags & LLT_CORRECT) {
		print_log("Triangulating and correcting %d frames...\n", *nframes);

#ifndef _OPENMP
		rvec *x2;
		snew(x2, 2 * *natoms);
#endif

#pragma omp parallel for shared(nframes,areas,x,natoms,flags)
		for(int i = 0; i < *nframes; ++i) {
#if defined _OPENMP && defined LLT_DEBUG
			print_log("%d threads triangulating.\n", omp_get_num_threads());
#endif
			real area1, area2, minx = FLT_MAX, maxx = FLT_MIN, dim;

			for(int j = 0; j < *natoms; ++j) {
				if(x[i][j][XX] < minx)	minx = x[i][j][XX];
				if(x[i][j][XX] > maxx)	maxx = x[i][j][XX];
			}
			dim = maxx - minx;

			area1 = tri_surface_area(x[i], *natoms, flags);

			// Correction for periodic bounds
			rvec *trans_x;
#ifdef _OPENMP
			snew(trans_x, 2 * *natoms); // if parallelized, each thread has its own array to write vectors
#else
			trans_x = x2; // otherwise, all iterations write to the same array
#endif
			memcpy(trans_x, x[i], sizeof(rvec) * *natoms);
			memcpy(trans_x + *natoms, x[i], sizeof(rvec) * *natoms);

			for(int j = *natoms; j < *natoms * 2; ++j) {
				trans_x[j][XX] += dim;
			}

			area2 = tri_surface_area(trans_x, 2 * *natoms, flags);

#ifdef _OPENMP
			sfree(trans_x);
#endif

			(*areas)[i] = 2 * area2 - 3 * area1;
		}
#ifndef _OPENMP
		sfree(x2);
#endif
	}
	else {
		print_log("Triangulating %d frames...\n", *nframes);

#pragma omp parallel for shared(nframes,areas,x,natoms,flags)
		for(int i = 0; i < *nframes; ++i) {
#if defined _OPENMP && defined LLT_DEBUG
			print_log("%d threads triangulating.\n", omp_get_num_threads());
#endif
			(*areas)[i] = tri_surface_area(x[i], *natoms, flags);
		}
	}

	// free memory
	for(int i = 0; i < *nframes; ++i) {
		sfree(x[i]);
	}
	sfree(x);
}


real tri_surface_area(rvec *x, int natoms, unsigned char flags) {
	static int iter = 0;

	char *tri_options = "zNQ";
	struct triangulateio tio;
	real area = 0;

	++iter;

	// input initialization
	snew(tio.pointlist, 2 * natoms);
	tio.pointmarkerlist = NULL;
	tio.numberofpoints = natoms;
	tio.numberofpointattributes = 0;

	for(int i = 0; i < natoms; ++i) {
		tio.pointlist[2*i] = x[i][0];
		tio.pointlist[2*i + 1] = x[i][1];
	}

	// output initialization
	tio.trianglelist = NULL;

	// triangulate the atoms
	triangulate(tri_options, &tio, &tio, NULL);

	if(flags & LLT_PRINT) { // print triangle data to files that can be viewed with triangle's 'showme' program
		char fname1[50], fname2[50];
		sprintf(fname1, "triangles%d.node", iter);
		sprintf(fname2, "triangles%d.ele", iter);
		print_trifiles(&tio, fname1, fname2);
	}

	sfree(tio.pointlist);

	// calculate surface area of triangles
	for(int i = 0; i < tio.numberoftriangles; ++i) {
		area += area_tri(x[tio.trianglelist[3*i]], x[tio.trianglelist[3*i + 1]], x[tio.trianglelist[3*i + 2]]);
	}

	sfree(tio.trianglelist);

	return area;
}


void print_areas(const char *fname, real *areas, int nframes, int natoms) {
	FILE *f = fopen(fname, "w");
	real sum = 0;

	fprintf(f, "FRAME\tAREA\tAREA/PARTICLE\n");
	for(int i = 0; i < nframes; ++i) {
		fprintf(f, "%d\t%f\t%f\n", i, areas[i], areas[i] / natoms);
		sum += areas[i];
	}

	fprintf(f, "\nAverage surface area: %f\n", sum / nframes);
	print_log("Average surface area: %f\n", sum / nframes);

	fprintf(f, "Average area per particle: %f\n", (sum / nframes) / natoms);
	print_log("Average area per particle: %f\n", (sum / nframes) / natoms);

	print_log("Surface areas saved to %s\n", fname);

	fclose(f);
}


void print_trifiles(struct triangulateio *tio, const char *node_name, const char *ele_name) {
	// print points to node file
	FILE *node = fopen(node_name, "w");

	fprintf(node, "%d\t2\t0\t0\n", tio->numberofpoints);
	for(int i = 0; i < tio->numberofpoints; ++i) {
		fprintf(node, "%d\t%f\t%f\n", i, tio->pointlist[2*i], tio->pointlist[2*i + 1]);
	}

	fclose(node);

	// print triangles to ele file
	FILE *ele = fopen(ele_name, "w");

	fprintf(ele, "%d\t3\t0\n", tio->numberoftriangles);
	for(int i = 0; i < tio->numberoftriangles; ++i) {
		fprintf(ele, "%d\t%d\t%d\t%d\n", 
			i, tio->trianglelist[3*i], tio->trianglelist[3*i + 1], tio->trianglelist[3*i + 2]);
	}

	fclose(ele);
}
