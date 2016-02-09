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
#include "gkut_io.h"
#include "gkut_log.h"
#include "smalloc.h"
#include "delaunay_tri.h"


void print_dtrifiles(	const struct dTriangulation *tri, 
						const char *node_name, 
						const char *ele_name);


void llt_delaunay_area(	const char *traj_fname, 
						const char *ndx_fname, 
						output_env_t *oenv, 
						real corr, 
						struct tri_area *areas, 
						unsigned char flags) {
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

#ifdef LLT_BENCH
	clock_t start = clock();
#endif

#ifdef _OPENMP
	if(!(flags & LLT_NOPAR))
		print_log("Triangulation will be parallelized.\n");
#endif

	// Calculate triangulated surface area for every frame
	dtinit(); // Initialize the delaunay triangulator
	snew(areas->area, areas->nframes);
	snew(areas->area2Dbox, areas->nframes);
	if(flags & LLT_2D)	snew(areas->area2D, areas->nframes);

	if(flags & LLT_CORRECT) { // add correction for periodic bounds
		print_log("Triangulating and correcting %d frames...\n", areas->nframes);

#pragma omp parallel for if(!(flags & LLT_NOPAR)) shared(areas,x,flags)
		for(int fr = 0; fr < areas->nframes; ++fr) {
			// Calculate number of edge points
			int n_edge_x = box[fr][0][0] / corr;
			int n_edge_y = box[fr][1][1] / corr;

			// Calculate box area
			areas->area2Dbox[fr] = box[fr][0][0] * box[fr][1][1];

			// z-coordinates of particles closest to box corners
			real bot_left = FLT_MAX, top_right = FLT_MIN, 
				top_left = FLT_MAX, bot_right = FLT_MIN, 
				avg_z;
			int bot_left_ind = 0, top_right_ind = 0, top_left_ind = 0, bot_right_ind = 0;
			
			// Find min max coordinates for each interval
			real *y_mins, *y_maxes, *x_mins, *x_maxes;
			snew(y_mins, n_edge_x + 1);
			snew(y_maxes, n_edge_x + 1);
			snew(x_mins, n_edge_y + 1);
			snew(x_maxes, n_edge_y + 1);

			int *y_min_inds, *y_max_inds, *x_min_inds, *x_max_inds;
			snew(y_min_inds, n_edge_x + 1);
			snew(y_max_inds, n_edge_x + 1);
			snew(x_min_inds, n_edge_y + 1);
			snew(x_max_inds, n_edge_y + 1);

			for(int i = 0; i <= n_edge_x; ++i)
				y_mins[i] = FLT_MAX;
			for(int i = 0; i <= n_edge_x; ++i)
				y_maxes[i] = FLT_MIN;
			for(int i = 0; i <= n_edge_y; ++i)
				x_mins[i] = FLT_MAX;
			for(int i = 0; i <= n_edge_y; ++i)
				x_maxes[i] = FLT_MIN;

			memset(y_min_inds, 0, sizeof(int) * (n_edge_x + 1));
			memset(y_max_inds, 0, sizeof(int) * (n_edge_x + 1));
			memset(x_min_inds, 0, sizeof(int) * (n_edge_y + 1));
			memset(x_max_inds, 0, sizeof(int) * (n_edge_y + 1));

			real dist, dY;
			int x_interval, y_interval;
			for(int j = 0; j < areas->natoms; ++j) {
				// min and max distance from origin
				dist = x[fr][j][XX] * x[fr][j][XX] + x[fr][j][YY] * x[fr][j][YY];
				if(dist < bot_left) {
					bot_left = dist;
					bot_left_ind = j;
				}
				if(dist > top_right) {
					top_right = dist;
					top_right_ind = j;
				}

				// min and max distance from top left corner
				dY = box[fr][1][1] - x[fr][j][YY];
				dist = x[fr][j][XX] * x[fr][j][XX] + dY * dY;
				if(dist < top_left) {
					top_left = dist;
					top_left_ind = j;
				}
				if(dist > bot_right) {
					bot_right = dist;
					bot_right_ind = j;
				}

				// Check min max y in x interval
				x_interval = (int)((x[fr][j][XX] / box[fr][0][0]) * n_edge_x);

				if(x[fr][j][YY] < y_mins[x_interval]) {
					y_mins[x_interval] = x[fr][j][YY];
					y_min_inds[x_interval] = j;
				}

				if(x[fr][j][YY] > y_maxes[x_interval]) {
					y_maxes[x_interval] = x[fr][j][YY];
					y_max_inds[x_interval] = j;
				}

				// Check min max x in y interval
				y_interval = (int)((x[fr][j][YY] / box[fr][1][1]) * n_edge_y);
				
				if(x[fr][j][XX] < x_mins[y_interval]) {
					x_mins[y_interval] = x[fr][j][XX];
					x_min_inds[y_interval] = j;
				}

				if(x[fr][j][XX] > x_maxes[y_interval]) {
					x_maxes[y_interval] = x[fr][j][XX];
					x_max_inds[y_interval] = j;
				}
			}

			sfree(y_mins);
			sfree(y_maxes);
			sfree(x_mins);
			sfree(x_maxes);

			avg_z = ( x[fr][bot_left_ind][ZZ] 
					+ x[fr][top_right_ind][ZZ] 
					+ x[fr][top_left_ind][ZZ] 
					+ x[fr][bot_right_ind][ZZ]) / 4.0;

			// add edge and corner points

			srenew(x[fr], areas->natoms + 2 * (n_edge_x + 1) + 2 * (n_edge_y + 1));
			int n = areas->natoms;

			// Add corner points
			x[fr][n][XX] 	= 0;
			x[fr][n][YY] 	= 0;
			x[fr][n++][ZZ] 	= avg_z;

			x[fr][n][XX] 	= box[fr][0][0];
			x[fr][n][YY] 	= 0;
			x[fr][n++][ZZ] 	= avg_z;

			x[fr][n][XX] 	= box[fr][0][0];
			x[fr][n][YY] 	= box[fr][1][1];
			x[fr][n++][ZZ] 	= avg_z;

			x[fr][n][XX] 	= 0;
			x[fr][n][YY] 	= box[fr][1][1];
			x[fr][n++][ZZ] 	= avg_z;

			// Add edge points
			real dist1, dist2;
			for(int j = 0; j < n_edge_x; ++j) {
				// Bottom edge
				x[fr][n][XX] = j * corr + corr / 2; // Go to middle of interval
				x[fr][n][YY] = 0;
				// edge Z coord is distance-from-edge-weighted average between the Zs of the two points closest to the two edges of this axis
				dist1 = x[fr][y_min_inds[j]][YY];
				dist2 = box[fr][1][1] - x[fr][y_max_inds[j]][YY];
				dist = dist1 + dist2;
				avg_z = x[fr][y_min_inds[j]][ZZ] - (dist1/dist)*(x[fr][y_min_inds[j]][ZZ]) 
					  + x[fr][y_max_inds[j]][ZZ] - (dist2/dist)*(x[fr][y_max_inds[j]][ZZ]);
				x[fr][n++][ZZ] = avg_z;

				// Top edge
				x[fr][n][XX] = j * corr + corr / 2;
				x[fr][n][YY] = box[fr][1][1];
				x[fr][n++][ZZ] = avg_z;
			}

			for(int j = 0; j < n_edge_y; ++j) {
				// Left edge
				x[fr][n][XX] = 0;
				x[fr][n][YY] = j * corr + corr / 2;
				
				dist1 = x[fr][x_min_inds[j]][XX];
				dist2 = box[fr][0][0] - x[fr][x_max_inds[j]][XX];
				dist = dist1 + dist2;
				avg_z = x[fr][x_min_inds[j]][ZZ] - (dist1/dist)*(x[fr][x_min_inds[j]][ZZ])
					  + x[fr][x_max_inds[j]][ZZ] - (dist2/dist)*(x[fr][x_max_inds[j]][ZZ]);
				x[fr][n++][ZZ] = avg_z;

				// Right edge
				x[fr][n][XX] = box[fr][0][0];
				x[fr][n][YY] = j * corr + corr / 2;
				x[fr][n++][ZZ] = avg_z;
			}

			sfree(y_min_inds);
			sfree(y_max_inds);
			sfree(x_min_inds);
			sfree(x_max_inds);

	#ifdef LLT_DEBUG
			FILE *f = fopen("points.txt", "w");

			for(int j = 0; j < n; ++j) {
				fprintf(f, "%d: %f\t%f\t%f\n", j, x[fr][j][XX], x[fr][j][YY], x[fr][j][ZZ]);
			}

			fclose(f);

			print_log("Points saved to points.txt for debugging.\n");
			exit(0);
	#endif

			// Calculate area including added edge and corner points
			real *a2D = NULL;
			if(flags & LLT_2D)	a2D = &(areas->area2D[fr]);
			delaunay_surface_area(x[fr], n, flags, a2D, &(areas->area[fr]));

			sfree(x[fr]);
		}
	}
	else { // triangulate without correction for periodic bounds
		print_log("Triangulating %d frames...\n", areas->nframes);

#pragma omp parallel for if(!(flags & LLT_NOPAR)) shared(areas,x,flags)
		for(int i = 0; i < areas->nframes; ++i) {
#if defined _OPENMP && defined LLT_DEBUG
			print_log("%d threads triangulating.\n", omp_get_num_threads());
#endif
			// 2D area of box
			areas->area2Dbox[i] = box[i][0][0] * box[i][1][1];

			real *a2D = NULL;
			if(flags & LLT_2D)	a2D = &(areas->area2D[i]);
			delaunay_surface_area(x[i], areas->natoms, flags, a2D, &(areas->area[i]));
			sfree(x[i]);
		}
	}

	sfree(x);
	sfree(box);

#ifdef LLT_BENCH
	clock_t clocks = clock() - start;
	print_log("Triangulation took %d clocks, %f seconds.\n", 
		clocks, (float)clocks/CLOCKS_PER_SEC);
#endif
}

void delaunay_surface_area(	const rvec *x, 
							int natoms, 
							unsigned char flags,
							real *a2D,
							real *a3D) {
	static int iter = 0;

	struct dTriangulation tri;
	++iter;

	// Input initialization
	snew(tri.points, 2 * natoms);
	tri.npoints = natoms;

	for(int i = 0; i < natoms; ++i) {
		tri.points[2*i] = x[i][XX];
		tri.points[2*i+1] = x[i][YY];
	}

	// triangulate
	dtriangulate(&tri);

	if(flags & LLT_PRINT) { // print triangle data to files that can be viewed with triangle's 'showme' program
		char fname1[50], fname2[50];
		sprintf(fname1, "triangles%d.node", iter);
		sprintf(fname2, "triangles%d.ele", iter);
		print_dtrifiles(&tri, fname1, fname2);
	} 

	sfree(tri.points);

	// calculate surface area of triangles
	if(a2D) {
		if(a3D) {
			*a2D = 0;
			*a3D = 0;
			rvec a, b, c;
			for(int i = 0; i < tri.ntriangles; ++i) {
				copy_rvec(x[tri.triangles[3*i]], a);
				copy_rvec(x[tri.triangles[3*i + 1]], b);
				copy_rvec(x[tri.triangles[3*i + 2]], c);

				(*a3D) += area_tri(a, b, c);

				a[ZZ] = 0;
				b[ZZ] = 0;
				c[ZZ] = 0;

				(*a2D) += area_tri(a, b, c);
			}
		}
		else {
			*a2D = 0;
			rvec a, b, c;
			for(int i = 0; i < tri.ntriangles; ++i) {
				copy_rvec(x[tri.triangles[3*i]], a);
				copy_rvec(x[tri.triangles[3*i + 1]], b);
				copy_rvec(x[tri.triangles[3*i + 2]], c);

				a[ZZ] = 0;
				b[ZZ] = 0;
				c[ZZ] = 0;

				(*a2D) += area_tri(a, b, c);
			}
		}
	}
	else if(a3D) {
		*a3D = 0;
		for(int i = 0; i < tri.ntriangles; ++i) {
			(*a3D) += area_tri( x[tri.triangles[3*i]], 
								x[tri.triangles[3*i + 1]], 
								x[tri.triangles[3*i + 2]]);
		}
	}

	free(tri.triangles);
}


void print_areas(const char *fname, const struct tri_area *areas) {
	FILE *f = fopen(fname, "w");
	real sum = 0;

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

	fclose(f);
	print_log("Surface areas saved to %s\n", fname);
}

void print_dtrifiles(	const struct dTriangulation *tri, 
						const char *node_name, 
						const char *ele_name) {
	// print points to node file
	FILE *node = fopen(node_name, "w");

	fprintf(node, "%d\t2\t0\t0\n", tri->npoints);
	for(int i = 0; i < tri->npoints; ++i) {
		fprintf(node, "%d\t%f\t%f\n", i, tri->points[2*i], tri->points[2*i + 1]);
	}

	fclose(node);

	// print triangles to ele file
	FILE *ele = fopen(ele_name, "w");

	fprintf(ele, "%d\t3\t0\n", tri->ntriangles);
	for(int i = 0; i < tri->ntriangles; ++i) {
		fprintf(ele, "%d\t%d\t%d\t%d\n", 
			i, tri->triangles[3*i], tri->triangles[3*i + 1], tri->triangles[3*i + 2]);
	}

	fclose(ele);
}

void free_tri_area(struct tri_area *areas) {
	if(areas->area)			sfree(areas->area);
	if(areas->area2D)		sfree(areas->area2D);
	if(areas->area2Dbox)	sfree(areas->area2Dbox);
}
