/*
 * Copyright Ahnaf Siddiqui and Sameer Varma.
 */

#include "ll_tessellation.h"

#include <float.h>
#include <stdio.h>
#include "smalloc.h"
#include "vec.h"

#include "gkut_io.h"
#include "gkut_log.h"

static real llt_diag, llt_diag2;

real weight_dist(rvec a, rvec b) {
	return llt_diag - sqrt(distance2(a, b));
}

real weight_dist2(rvec a, rvec b) {
	return llt_diag2 - distance2(a, b);
}


void llt_area(const char *traj_fname, const char *ndx_fname, 
	real cell_width, real (*fweight)(rvec, rvec), output_env_t *oenv, struct tessellated_grid *grid) {
	rvec **pre_x, **x;
	int nframes, natoms;

	read_traj(traj_fname, &pre_x, &nframes, &natoms, oenv);

	// Filter trajectory by index file if present
	if(ndx_fname != NULL) {
		const int NUMGROUPS = 1;
		int *isize;
		atom_id **indx;
		char **grp_names;

		snew(isize, NUMGROUPS);
		snew(indx, NUMGROUPS);
		snew(grp_names, NUMGROUPS);

		rd_index(ndx_fname, NUMGROUPS, isize, indx, grp_names);
		sfree(grp_names);

		natoms = isize[0];
		sfree(isize);

		snew(x, nframes);
		for(int i = 0; i < nframes; ++i) {
			snew(x[i], natoms);
			for(int j = 0; j < natoms; ++j) {
				copy_rvec(pre_x[i][indx[0][j]], x[i][j]);
			}
		}

		// free memory
		sfree(indx[0]);
		sfree(indx);

		for(int i = 0; i < nframes; ++i) {
			sfree(pre_x[i]);
		}
		sfree(pre_x);
	}
	else {
		x = pre_x;
	}

	f_llt_area(x, nframes, natoms, cell_width, fweight, grid);

	// free memory
	for(int i = 0; i < nframes; ++i) {
		sfree(x[i]);
	}
	sfree(x);
}


void f_llt_area(rvec **x, int nframes, int natoms, 
	real cell_width, real (*fweight)(rvec, rvec), struct tessellated_grid *grid) {
	construct_grid(x, nframes, natoms, cell_width, grid);

	load_grid(x, nframes, natoms, fweight, grid);

	gen_heightmap(grid);

	tessellate_area(grid);

	grid->area_per_particle = grid->surface_area / natoms;
}


void construct_grid(rvec **x, int nframes, int natoms, real cell_width, struct tessellated_grid *grid) {
	real minx = FLT_MAX, miny = FLT_MAX, minz = FLT_MAX, 
		maxx = FLT_MIN, maxy = FLT_MIN, maxz = FLT_MIN;
	int dimx, dimy, dimz;

	for(int fr = 0; fr < nframes; ++fr) {
		for(int a = 0; a < natoms; ++a) {
			if(x[fr][a][XX] < minx)	minx = x[fr][a][XX];
			if(x[fr][a][XX] > maxx)	maxx = x[fr][a][XX];

			if(x[fr][a][YY] < miny)	miny = x[fr][a][YY];
			if(x[fr][a][YY] > maxy)	maxy = x[fr][a][YY];

			if(x[fr][a][ZZ] < minz)	minz = x[fr][a][ZZ];
			if(x[fr][a][ZZ] > maxz)	maxz = x[fr][a][ZZ];
		}
	}

	// # weights in each dim is the # grid cells - 1 + an extra grid cell (bc of int cast floor) + 1 for the last grid point
	dimx = ((int)((maxx - minx)/cell_width) + 2);
	dimy = ((int)((maxy - miny)/cell_width) + 2);
	dimz = ((int)((maxz - minz)/cell_width) + 2);

	snew(grid->weights, dimx * dimy * dimz);
	snew(grid->heightmap, dimx * dimy);
	snew(grid->areas, (dimx-1) * (dimy-1));
	grid->dimx = dimx, grid->dimy = dimy, grid->dimz = dimz;
	grid->cell_width = cell_width;
	grid->minx = minx, grid->miny = miny, grid->minz = minz;
#ifdef LLT_DEBUG
	print_log("maxx = %f, maxy = %f, maxz = %f\n", maxx, maxy, maxz);
#endif
}


void load_grid(rvec **x, int nframes, int natoms, real (*fweight)(rvec, rvec), struct tessellated_grid *grid) {
	real *weights = grid->weights;
	int dimx = grid->dimx, dimy = grid->dimy, dimz = grid->dimz;
	int dimyz = dimy * dimz;
	real cell_width = grid->cell_width;
	real minx = grid->minx, miny = grid->miny, minz = grid->minz;

	llt_diag2 = 3 * cell_width * cell_width;
	llt_diag = sqrt(llt_diag2);

	rvec grid_point;
	int xi, yi, zi;

	for(int fr = 0; fr < nframes; ++fr) {
		for(int a = 0; a < natoms; ++a) {
			// Indices of the origin point of the grid cell surrounding this atom
			xi = (int)((x[fr][a][XX] - minx)/cell_width);
			yi = (int)((x[fr][a][YY] - miny)/cell_width);
			zi = (int)((x[fr][a][ZZ] - minz)/cell_width);

			// Load the eight grid points around this atom. Closer distance to atom = higher weight
			grid_point[XX] = minx + xi * cell_width;
			grid_point[YY] = miny + yi * cell_width; 
			grid_point[ZZ] = minz + zi * cell_width;
			// This order of operations is an attempt to minimize cache misses
			weights[xi*dimyz + yi*dimz + zi] 			+= fweight(x[fr][a], grid_point);
			grid_point[ZZ] += cell_width;
			weights[xi*dimyz + yi*dimz + zi+1] 			+= fweight(x[fr][a], grid_point);
			grid_point[YY] += cell_width;
			grid_point[ZZ] -= cell_width;
			weights[xi*dimyz + (yi+1)*dimz + zi] 		+= fweight(x[fr][a], grid_point);
			grid_point[ZZ] += cell_width;
			weights[xi*dimyz + (yi+1)*dimz + zi+1] 		+= fweight(x[fr][a], grid_point);
			grid_point[XX] += cell_width;
			grid_point[YY] -= cell_width;
			grid_point[ZZ] -= cell_width;
			weights[(xi+1)*dimyz + yi*dimz + zi] 		+= fweight(x[fr][a], grid_point);
			grid_point[ZZ] += cell_width;
			weights[(xi+1)*dimyz + yi*dimz + zi+1] 		+= fweight(x[fr][a], grid_point);
			grid_point[YY] += cell_width;
			grid_point[ZZ] -= cell_width;
			weights[(xi+1)*dimyz + (yi+1)*dimz + zi] 	+= fweight(x[fr][a], grid_point);
			grid_point[ZZ] += cell_width;
			weights[(xi+1)*dimyz + (yi+1)*dimz + zi+1] 	+= fweight(x[fr][a], grid_point);
		}
	}
}


void gen_heightmap(struct tessellated_grid *grid) {
	real *weights = grid->weights;
	int dimx = grid->dimx, dimy = grid->dimy, dimz = grid->dimz;
	int dimyz = dimy * dimz;

	int *heightmap = grid->heightmap;

	int maxz;
	real max_weight;
	int num_empty = 0;

	for(int x = 0; x < dimx; ++x) {
		for(int y = 0; y < dimy; ++y) {
			maxz = -1; // If none of the weights in this column is > ~0, then this column's z index will be -1 
			max_weight = 2 * FLT_EPSILON; // To protect the criterion above in case of floating point imprecision
			for(int z = 0; z < dimz; ++z) {
				if(weights[x*dimyz + y*dimz + z] > max_weight) {
					max_weight = weights[x*dimyz + y*dimz + z];
					maxz = z;
				}
			}
			heightmap[x*dimy + y] = maxz;
			num_empty += maxz < 0;
		}
	}

	grid->num_empty = num_empty;
}


void tessellate_area(struct tessellated_grid *grid) {
	int dimx = grid->dimx, dimy = grid->dimy, dimz = grid->dimz;
	int *heightmap = grid->heightmap;
	real *areas = grid->areas;
	real cell_width = grid->cell_width;

	real tot_area = 0, cell_area = 0;

	int i_heights[4];
	rvec corners[] = {
		{0.0, 0.0, 0.0}, 
		{0.0, cell_width, 0.0}, 
		{cell_width, 0.0, 0.0}, 
		{cell_width, cell_width, 0.0}
	};
	rvec ab, ac, ad, cpr;

#ifdef LLT_DEBUG
			print_log("Corner height indices:\n");
#endif

	for(int x = 0; x < dimx - 1; ++x) {
		for(int y = 0; y < dimy - 1; ++y) {
			if(((i_heights[0] = heightmap[x*dimy + y]) < 0) 
				|| ((i_heights[1] = heightmap[x*dimy + y+1]) < 0) 
				|| ((i_heights[2] = heightmap[(x+1)*dimy + y]) < 0) 
				|| ((i_heights[3] = heightmap[(x+1)*dimy + y+1]) < 0)) {
				areas[x*dimy + y] = 0.0;
				continue;
			}

#ifdef LLT_DEBUG
			print_log("Cell [%d][%d]: %d %d %d %d\n", 
				x, y, i_heights[0], i_heights[1], i_heights[2], i_heights[3]);
#endif
			corners[0][ZZ] = i_heights[0] * cell_width;
			corners[1][ZZ] = i_heights[1] * cell_width;
			corners[2][ZZ] = i_heights[2] * cell_width;
			corners[3][ZZ] = i_heights[3] * cell_width;

			rvec_sub(corners[1], corners[0], ab);
			rvec_sub(corners[2], corners[0], ac);
			rvec_sub(corners[3], corners[0], ad);

			cprod(ab, ad, cpr);
			cell_area = norm(cpr) / 2.0;

			cprod(ad, ac, cpr);
			cell_area += norm(cpr) / 2.0;

			tot_area += cell_area;

			areas[x*(dimy-1) + y] = cell_area;
		}
	}

	grid->surface_area = tot_area;
}


void print_grid(struct tessellated_grid *grid, const char *fname) {
	int dimx = grid->dimx, dimy = grid->dimy, dimz = grid->dimz;

	FILE *f = fopen(fname, "w");

	fprintf(f, "Grid points: dimx = %d, dimy = %d, dimz = %d\n", dimx, dimy, dimz);
	fprintf(f, "Grid cell width = %f\n", grid->cell_width);
	fprintf(f, "Trajectory origin: minx = %f, miny = %f, minz = %f\n", grid->minx, grid->miny, grid->minz);

	fprintf(f, "\nWeights ([x][y]: z weights):");
	for(int x = 0; x < dimx; ++x) {
		for(int y = 0; y < dimy; ++y) {
			fprintf(f, "\n[%d][%d]: ", x, y);
			for(int z = 0; z < dimz; ++z) {
				fprintf(f, "%f ", grid->weights[x*dimy*dimz + y*dimz + z]);
			}
		}
	}

	fprintf(f, "\n\nHeightmap ([x]: y heights in z index. -1 means empty.):");
	for(int x = 0; x < dimx; ++x) {
		fprintf(f, "\n[%d]: ", x);
		for(int y = 0; y < dimy; ++y) {
			fprintf(f, "%d ", grid->heightmap[x*dimy + y]);
		}
	}

	fprintf(f, "\n\nGrid cell areas (Cell with origin indices [x][y]: area):\n");
	for(int x = 0; x < dimx - 1; ++x) {
		for(int y = 0; y < dimy - 1; ++y) {
			fprintf(f, "Cell [%d][%d]: %f\n", x, y, grid->areas[x*(dimy-1) + y]);
		}
	}

	fprintf(f, "\n%d grid cell(s) have empty (zero-weight) corner(s) and are excluded from tessellation.\n", grid->num_empty);
	fprintf(f, "Total tessellated surface area: %f\n", grid->surface_area);
	fprintf(f, "Tessellated surface area per particle: %f\n", grid->area_per_particle);

	fclose(f);
}


void free_grid(struct tessellated_grid *grid) {
	sfree(grid->weights);
	sfree(grid->heightmap);
	sfree(grid->areas);
}
