/*
 * Copyright Ahnaf Siddiqui and Sameer Varma.
 */

#include "llt_tri.h"

#include "gkut_io.h"
#include "gkut_log.h"
#include "smalloc.h"


void llt_tri_area(const char *traj_fname, const char *ndx_fname, output_env_t *oenv, 
	real **areas, int *nframes, int *natoms) {
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
	for(int i = 0; i < *nframes; ++i) {
		(*areas)[i] = tri_surface_area(x[i], *natoms);
	}

	// free memory
	for(int i = 0; i < *nframes; ++i) {
		sfree(x[i]);
	}
	sfree(x);
}


real tri_surface_area(rvec *x, int natoms) {
#ifdef LLT_DEBUG
	static int iter = 0;
	++iter;
#endif
	char *tri_options = "zNQ";
	struct triangulateio tio;
	real area = 0;

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

#ifdef LLT_DEBUG
	char fname1[50], fname2[50];
	sprintf(fname1, "triangles%d.node", iter);
	sprintf(fname2, "triangles%d.ele", iter);
	print_trifiles(&tio, fname1, fname2);
#endif

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
