/*
 * Copyright Ahnaf Siddiqui and Sameer Varma.
 */

#include "ll_tessellation.h"

#include <float.h>
#include <stdio.h>
#include "smalloc.h"
#include "vec.h"

#include "gkut_io.h"

// #define FRAMESTEP 500 // The number of new frames by which to reallocate an array of length # trajectory frames

// static void read_traj(const char *traj_fname, rvec ***x, int *nframes, int *natoms, output_env_t *oenv) {
// 	t_trxstatus *status = NULL;
// 	real t;
// 	matrix box;
// 	int est_frames = FRAMESTEP;
// 	*nframes = 0;

// 	snew(*x, est_frames);
// 	*natoms = read_first_x(*oenv, &status, traj_fname, &t, &((*x)[0]), box);

// 	do {
// 		(*nframes)++;
// 		if(*nframes >= est_frames) {
// 			est_frames += FRAMESTEP;
// 			srenew(*x, est_frames);
// 		}
// 		snew((*x)[*nframes], *natoms);
// 	} while(read_next_x(*oenv, status, &t,
// #ifndef GRO_V5 
// 		*natoms,
// #endif
// 		(*x)[*nframes], box));

// 	sfree((*x)[*nframes]); // Nothing was read to the last allocated frame
// 	close_trx(status);
// }


real tessellate_area(const char *traj_fname, const char *ndx_fname, int numcells, output_env_t *oenv) {
	rvec **pre_x, **x;
	int nframes, natoms;
	struct weighted_grid grid;

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

	construct_grid(x, nframes, natoms, numcells, &grid);

#ifdef LLT_DEBUG
	printf("Grid: \n");
	printf("dimx = %d, dimy = %d, dimz = %d\n", grid.dimx, grid.dimy, grid.dimz);
	printf("cell width = %f\n", grid.cell_width);
	printf("minx = %f, miny = %f, minz = %f\n", grid.minx, grid.miny, grid.minz);
#endif

	load_grid(x, nframes, natoms, &grid);

	// free memory
	free_grid(&grid);

	for(int i = 0; i < nframes; ++i) {
		sfree(x[i]);
	}
	sfree(x);

	return 0;
}

void construct_grid(rvec **x, int nframes, int natoms, int numcells, struct weighted_grid *grid) {
	real minx = FLT_MAX, miny = FLT_MAX, minz = FLT_MAX, 
		maxx = FLT_MIN, maxy = FLT_MIN, maxz = FLT_MIN;
	real distx, disty, distz, cell_width;
	int dimx, dimz, dimy;

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

	distx = maxx - minx, disty = maxy - miny, distz = maxz - minz;
	real distb;
	cell_width = (distx > (distb = disty > distz ? disty : distz) ? distx : distb) / numcells;

	// # weights in each dim is the # grid cells - 1 + an extra grid cell (bc of int cast floor) + 1 for the last grid point
	dimx = ((int)(distx/cell_width) + 2), dimz = ((int)(distz/cell_width) + 2), dimy = ((int)(disty/cell_width) + 2);

	snew(grid->weights, dimx * dimy * dimz);
	grid->dimx = dimx, grid->dimy = dimy, grid->dimz = dimz;
	grid->cell_width = cell_width;
	grid->minx = minx, grid->miny = miny, grid->minz = minz;
#ifdef LLT_DEBUG
	printf("maxx = %f, maxy = %f, maxz = %f\n", maxx, maxy, maxz);
#endif
}

void load_grid(rvec **x, int nframes, int natoms, struct weighted_grid *grid) {
	real *weights = grid->weights;
	int dimx = grid->dimx, dimy = grid->dimy, dimz = grid->dimz;
	real cell_width = grid->cell_width;
	real minx = grid->minx, miny = grid->miny, minz = grid->minz;

	real diag_sq = 3 * cell_width * cell_width;
	rvec grid_point;
	int xi, yi, zi;
	int dimyz = dimy * dimz;

	for(int fr = 0; fr < nframes; ++fr) {
		for(int a = 0; a < natoms; ++a) {
			// Indices of the origin point of the grid cell surrounding this atom
			xi = (int)((x[fr][a][XX] - minx)/cell_width);
			yi = (int)((x[fr][a][YY] - miny)/cell_width);
			zi = (int)((x[fr][a][ZZ] - minz)/cell_width);

			// Load the eight grid points around this atom. Closer distance to atom = higher weight
			grid_point[XX] = minx + xi * cell_width, 
				grid_point[YY] = miny + yi * cell_width, 
				grid_point[ZZ] = minz + zi * cell_width;
			*(weights + xi * dimyz + yi * dimz + zi) += diag_sq - distance2(x[fr][a], grid_point);
			grid_point[XX] += cell_width;
			*(weights + (xi+1) * dimyz + yi * dimz + zi) += diag_sq - distance2(x[fr][a], grid_point);
			grid_point[YY] += cell_width;
			*(weights + (xi+1) * dimyz + (yi+1) * dimz + zi) += diag_sq - distance2(x[fr][a], grid_point);
			grid_point[XX] -= cell_width;
			*(weights + xi * dimyz + (yi+1) * dimz + zi) += diag_sq - distance2(x[fr][a], grid_point);
			grid_point[ZZ] += cell_width;
			*(weights + xi * dimyz + (yi+1) * dimz + zi+1) += diag_sq - distance2(x[fr][a], grid_point);
			grid_point[YY] -= cell_width;
			*(weights + xi * dimyz + yi * dimz + zi+1) += diag_sq - distance2(x[fr][a], grid_point);
			grid_point[XX] += cell_width;
			*(weights + (xi+1) * dimyz + yi * dimz + zi+1) += diag_sq - distance2(x[fr][a], grid_point);
			grid_point[YY] += cell_width;
			*(weights + (xi+1) * dimyz + (yi+1) * dimz + zi+1) += diag_sq - distance2(x[fr][a], grid_point);
		}
	}

#ifdef LLT_DEBUG
	printf("Weights: \n");
	for(int i = 0; i < dimx * dimyz; i++) {
		printf("%f ", *(weights + i));
	}
	printf("\n");
#endif
}

void free_grid(struct weighted_grid *grid) {
	sfree(grid->weights);
}
