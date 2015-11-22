/*
 * Copyright Ahnaf Siddiqui and Sameer Varma.
 */

#include "trt_tessellation.h"

#include <float.h>
#include <stdio.h>
#include "smalloc.h"
#include "vec.h"

#define FRAMESTEP 500 // The number of new frames by which to reallocate an array of length # trajectory frames


static void read_traj(const char *traj_fname, rvec ***x, int *nframes, int *natoms, output_env_t *oenv) {
	t_trxstatus *status = NULL;
	real t;
	matrix box;
	int est_frames = FRAMESTEP;
	*nframes = 0;

	snew(*x, est_frames);
	*natoms = read_first_x(*oenv, &status, traj_fname, &t, &((*x)[0]), box);

	do {
		(*nframes)++;
		if(*nframes >= est_frames) {
			est_frames += FRAMESTEP;
			srenew(*x, est_frames);
		}
		snew((*x)[*nframes], *natoms);
	} while(read_next_x(*oenv, status, &t,
#ifndef GRO_V5 
		*natoms,
#endif
		(*x)[*nframes], box));

	sfree((*x)[*nframes]); // Nothing was read to the last allocated frame
	close_trx(status);
}


real tr_tessellate_area(const char *traj_fname, const char *ndx_fname, int numcells, output_env_t *oenv) {
	rvec **x;
	int nframes, natoms;
	struct weighted_grid grid;

	read_traj(traj_fname, &x, &nframes, &natoms, oenv);

	construct_grid(x, nframes, natoms, numcells, &grid);

	// printf("dimx = %d, dimz = %d, dimy = %d\n", grid.dimx, grid.dimz, grid.dimy);
	// printf("cell width = %f\n", grid.cell_width);
	// printf("minx = %f, minz = %f, miny = %f\n", grid.minx, grid.minz, grid.miny);

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

	rvec cur_x;
	for(int fr = 0; fr < nframes; ++fr) {
		for(int a = 0; a < natoms; ++a) {
			copy_rvec(x[fr][a], cur_x);
			if(cur_x[XX] < minx)	minx = cur_x[XX];
			if(cur_x[XX] > maxx)	maxx = cur_x[XX];

			if(cur_x[YY] < miny)	miny = cur_x[YY];
			if(cur_x[YY] > maxy)	maxy = cur_x[YY];

			if(cur_x[ZZ] < minz)	minz = cur_x[ZZ];
			if(cur_x[ZZ] > maxz)	maxz = cur_x[ZZ];
		}
	}

	distx = maxx - minx, disty = maxy - miny, distz = maxz - minz;
	real distb;
	cell_width = (distx > (distb = disty > distz ? disty : distz) ? distx : distb) / numcells;

	// # weights in each dim is the # grid cells - 1 + an extra grid cell (bc of int cast floor) + 1 for the last grid point
	dimx = ((int)(distx/cell_width) + 2), dimz = ((int)(distz/cell_width) + 2), dimy = ((int)(disty/cell_width) + 2);

	snew(grid->weights, dimx * dimz * dimy);
	grid->dimx = dimx, grid->dimz = dimz, grid->dimy = dimy;
	grid->cell_width = cell_width;
	grid->minx = minx, grid->minz = minz, grid->miny = miny;
	// printf("maxx = %f, maxz = %f, maxy = %f\n", maxx, maxz, maxy);
}

void load_grid(rvec **x, int nframes, int natoms, struct weighted_grid *grid) {
	real *weights = grid->weights;
	int dimx = grid->dimx, dimz = grid->dimz, dimy = grid->dimy;
	real cell_width = grid->cell_width;
	real minx = grid->minx, miny = grid->miny, minz = grid->minz;

	real diag_sq = 3 * cell_width * cell_width;
	rvec cur_x, grid_point;
	int xi, yi, zi;

	for(int fr = 0; fr < nframes; ++fr) {
		for(int a = 0; a < natoms; ++a) {
			copy_rvec(x[fr][a], cur_x);

			// Indices of the origin point of the grid cell surrounding this atom
			xi = (int)((cur_x[XX] - minx)/cell_width);
			yi = (int)((cur_x[YY] - miny)/cell_width);
			zi = (int)((cur_x[ZZ] - minz)/cell_width);

			// Load the eight grid points around this atom. Closer distance to atom = higher weight
			grid_point[XX] = minx + xi * cell_width, 
				grid_point[YY] = miny + yi * cell_width, 
				grid_point[ZZ] = minz + zi * cell_width;
			*(weights + xi * dimz * dimy + zi * dimy + yi) += diag_sq - distance2(cur_x, grid_point);
			grid_point[XX] += cell_width;
			*(weights + (xi+1) * dimz * dimy + zi * dimy + yi) += diag_sq - distance2(cur_x, grid_point);
			grid_point[ZZ] += cell_width;
			*(weights + (xi+1) * dimz * dimy + (zi+1) * dimy + yi) += diag_sq - distance2(cur_x, grid_point);
			grid_point[XX] -= cell_width;
			*(weights + xi * dimz * dimy + (zi+1) * dimy + yi) += diag_sq - distance2(cur_x, grid_point);
			grid_point[YY] += cell_width;
			*(weights + xi * dimz * dimy + (zi+1) * dimy + yi+1) += diag_sq - distance2(cur_x, grid_point);
			grid_point[ZZ] -= cell_width;
			*(weights + xi * dimz * dimy + zi * dimy + yi+1) += diag_sq - distance2(cur_x, grid_point);
			grid_point[XX] += cell_width;
			*(weights + (xi+1) * dimz * dimy + zi * dimy + yi+1) += diag_sq - distance2(cur_x, grid_point);
			grid_point[ZZ] += cell_width;
			*(weights + (xi+1) * dimz * dimy + (zi+1) * dimy + yi+1) += diag_sq - distance2(cur_x, grid_point);
		}
	}

	for(int i = 0; i < dimx * dimz * dimy; i++) {
		printf("%f ", *(weights + i));
	}
	printf("\n");
}

void free_grid(struct weighted_grid *grid) {
	// for(int i = 0; i < grid->dimx; ++i) {
	// 	for(int j = 0; j < grid->dimz; ++j) {
	// 		sfree(grid->weights[i][j]);
	// 	}
	// 	sfree(grid->weights[i]);
	// }
	// sfree(grid->weights);
	sfree(grid->weights);
}
