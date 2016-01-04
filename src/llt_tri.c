/*
 * Copyright Ahnaf Siddiqui and Sameer Varma.
 */

#include "llt_tri.h"

#include <float.h>
#include <string.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#ifdef LLT_BENCH
#include <time.h>
#endif
#include "triangle.h"
#include "gkut_io.h"
#include "gkut_log.h"
#include "smalloc.h"

void print_trifiles(struct triangulateio *tio, const char *node_name, const char *ele_name);


void llt_tri_area(const char *traj_fname, const char *ndx_fname, output_env_t *oenv, 
	struct tri_area *areas, unsigned char flags) {
	rvec **pre_x, **x;
	matrix *box;

	areas->area = NULL;
	areas->area2D = NULL;
	areas->area2Dbox = NULL;
	areas->area1 = NULL;
	areas->area2 = NULL;

	read_traj(traj_fname, &pre_x, &box, &(areas->nframes), &(areas->natoms), oenv);

	// Filter trajectory by index file if present
	if(ndx_fname != NULL) {
		ndx_filter_traj(ndx_fname, pre_x, &x, areas->nframes, &(areas->natoms));

		for(int i = 0; i < areas->nframes; ++i) {
			sfree(pre_x[i]);
		}
		sfree(pre_x);
	}
	else {
		x = pre_x;
	}

#ifdef LLT_DEBUG
	// test filtering
	print_traj(x, areas->nframes, areas->natoms, "traj.dat");
#endif

#ifdef LLT_BENCH
	clock_t start = clock();
#endif

#ifdef _OPENMP
	print_log("Triangulation will be parallelized.\n");
#endif
	// Calculate triangulated surface area for every frame
	snew(areas->area, areas->nframes);
	snew(areas->area2Dbox, areas->nframes);
	if(flags & LLT_2D)	snew(areas->area2D, areas->nframes);

	if(flags & LLT_CORRECT) {
		print_log("Triangulating and correcting %d frames...\n", areas->nframes);
		snew(areas->area1, areas->nframes);
		snew(areas->area2, areas->nframes);

#ifndef _OPENMP
		rvec *x2;
		snew(x2, 2 * areas->natoms);
#endif

#pragma omp parallel for shared(areas,x,flags)
		for(int i = 0; i < areas->nframes; ++i) {
#if defined _OPENMP && defined LLT_DEBUG
			print_log("%d threads triangulating.\n", omp_get_num_threads());
#endif
			
			// 2D area of box
			areas->area2Dbox[i] = box[i][0][0] * box[i][1][1];

			real *a2D = NULL;
			if(flags & LLT_2D)	a2D = &(areas->area2D[i]);
			// Get area without correction
			areas->area1[i] = tri_surface_area(x[i], areas->natoms, flags, a2D);

			// Correction for periodic bounds
			rvec *trans_x;
#ifdef _OPENMP
			snew(trans_x, 2 * areas->natoms); // if parallelized, each thread has its own array to write vectors
#else
			trans_x = x2; // otherwise, all iterations write to the same array
#endif
			// Generate mirror image of particles...
			memcpy(trans_x, x[i], sizeof(rvec) * areas->natoms);
			memcpy(trans_x + areas->natoms, x[i], sizeof(rvec) * areas->natoms);
			// ...and translate it by box width
			for(int j = areas->natoms; j < areas->natoms * 2; ++j) {
				trans_x[j][XX] += box[i][0][0];
			}

			// Get area with translated image
			areas->area2[i] = tri_surface_area(trans_x, 2 * areas->natoms, flags, NULL);

#ifdef _OPENMP
			sfree(trans_x);
#endif

			// Corrected area = A1 + 2*(A2 - 2*A1)
			areas->area[i] = 2 * areas->area2[i] - 3 * areas->area1[i];
			// areas->area[i] = areas->area1[i] + 2 * (areas->area2[i] - 2 * areas->area1[i]);
		}
#ifndef _OPENMP
		sfree(x2);
#endif
	}
	else {
		print_log("Triangulating %d frames...\n", areas->nframes);

#pragma omp parallel for shared(areas,x,flags)
		for(int i = 0; i < areas->nframes; ++i) {
#if defined _OPENMP && defined LLT_DEBUG
			print_log("%d threads triangulating.\n", omp_get_num_threads());
#endif
			// 2D area of box
			areas->area2Dbox[i] = box[i][0][0] * box[i][1][1];

			real *a2D = NULL;
			if(flags & LLT_2D)	a2D = &(areas->area2D[i]);
			areas->area[i] = tri_surface_area(x[i], areas->natoms, flags, a2D);
		}
	}

	// free memory
	for(int i = 0; i < areas->nframes; ++i) {
		sfree(x[i]);
	}
	sfree(x);
	sfree(box);

#ifdef LLT_BENCH
	clock_t clocks = clock() - start;
	print_log("Triangulation took %d clocks, %f seconds.\n", clocks, (float)clocks/CLOCKS_PER_SEC);
#endif
}


real tri_surface_area(rvec *x, int natoms, unsigned char flags, real *a2D) {
	static int iter = 0;

	char *tri_options = "zNQ";
	struct triangulateio tio;
	real a3D = 0;

	++iter;

	// input initialization
	snew(tio.pointlist, 2 * natoms);
	tio.pointmarkerlist = NULL;
	tio.numberofpoints = natoms;
	tio.numberofpointattributes = 0;

	for(int i = 0; i < natoms; ++i) {
		tio.pointlist[2*i] = x[i][XX];
		tio.pointlist[2*i + 1] = x[i][YY];
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
	if(a2D != NULL) {
		*a2D = 0;
		rvec a, b, c;
		for(int i = 0; i < tio.numberoftriangles; ++i) {
			copy_rvec(x[tio.trianglelist[3*i]], a);
			copy_rvec(x[tio.trianglelist[3*i + 1]], b);
			copy_rvec(x[tio.trianglelist[3*i + 2]], c);

			a3D += area_tri(a, b, c);

			a[ZZ] = 0;
			b[ZZ] = 0;
			c[ZZ] = 0;

			(*a2D) += area_tri(a, b, c);
		}
	}
	else {
		for(int i = 0; i < tio.numberoftriangles; ++i) {
			a3D += area_tri(x[tio.trianglelist[3*i]], x[tio.trianglelist[3*i + 1]], x[tio.trianglelist[3*i + 2]]);
		}
	}

	sfree(tio.trianglelist);

	return a3D;
}

// TODO: make more flexible printing conditionals
void print_areas(const char *fname, struct tri_area *areas) {
	FILE *f = fopen(fname, "w");
	real sum = 0;

	if(areas->area1 && areas->area2) { // Corrected for periodic bounding conditions
		if(areas->area2D) {
			fprintf(f, "# FRAME\tAREA1\tAREA2\tCORRECTED-AREA\t2DAREA\tBOX-AREA\t\"\"/particle\n");
			for(int i =0; i < areas->nframes; ++i) {
				fprintf(f, "%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", i, 
					areas->area1[i], areas->area2[i], areas->area[i], areas->area2D[i], areas->area2Dbox[i],
					areas->area1[i] / areas->natoms, areas->area2[i] / areas->natoms, 
					areas->area[i] / areas->natoms, areas->area2D[i] / areas->natoms,
					areas->area2Dbox[i] / areas->natoms);
				sum += areas->area[i];
			}
		}
		else {
			fprintf(f, "# FRAME\tAREA1\tAREA2\tCORRECTED-AREA\tBOX-AREA\t\"\"/particle\n");
			for(int i =0; i < areas->nframes; ++i) {
				fprintf(f, "%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", i, 
					areas->area1[i], areas->area2[i], areas->area[i], areas->area2Dbox[i],
					areas->area1[i] / areas->natoms, areas->area2[i] / areas->natoms, 
					areas->area[i] / areas->natoms, areas->area2Dbox[i] / areas->natoms);
				sum += areas->area[i];
			}
		}

		fprintf(f, "\n# Average corrected surface area: %f\n", sum / areas->nframes);
		print_log("Average corrected surface area: %f\n", sum / areas->nframes);

		fprintf(f, "# Average corrected area per particle: %f\n", (sum / areas->nframes) / areas->natoms);
		print_log("Average corrected area per particle: %f\n", (sum / areas->nframes) / areas->natoms);
	}
	else { // Not corrected
		if(areas->area2D) {
			fprintf(f, "# FRAME\tAREA\t2DAREA\tBOX-AREA\t\"\"/PARTICLE\n");
			for(int i = 0; i < areas->nframes; ++i) {
				fprintf(f, "%d\t%f\t%f\t%f\t%f\t%f\t%f\n", i, areas->area[i], areas->area2D[i], areas->area2Dbox[i], 
					areas->area[i] / areas->natoms, areas->area2D[i] / areas->natoms, areas->area2Dbox[i] / areas->natoms);
				sum += areas->area[i];
			}
		}
		else {
			fprintf(f, "# FRAME\tAREA\tBOX-AREA\t\"\"/PARTICLE\n");
			for(int i = 0; i < areas->nframes; ++i) {
				fprintf(f, "%d\t%f\t%f\t%f\t%f\n", i, areas->area[i], areas->area2Dbox[i], 
					areas->area[i] / areas->natoms, areas->area2Dbox[i] / areas->natoms);
				sum += areas->area[i];
			}
		}

		fprintf(f, "\n# Average surface area: %f\n", sum / areas->nframes);
		print_log("Average surface area: %f\n", sum / areas->nframes);

		fprintf(f, "# Average area per particle: %f\n", (sum / areas->nframes) / areas->natoms);
		print_log("Average area per particle: %f\n", (sum / areas->nframes) / areas->natoms);
	}

	fclose(f);
	print_log("Surface areas saved to %s\n", fname);
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

void free_tri_area(struct tri_area *areas) {
	if(areas->area)			sfree(areas->area);
	if(areas->area2D)		sfree(areas->area2D);
	if(areas->area2Dbox)	sfree(areas->area2Dbox);
	if(areas->area1)		sfree(areas->area1);
	if(areas->area2)		sfree(areas->area2);
}
