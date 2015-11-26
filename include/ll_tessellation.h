/*
 * Copyright Ahnaf Siddiqui and Sameer Varma.
 */

#ifndef LL_TESSELLATION_H
#define LL_TESSELLATION_H

#ifdef GRO_V5
#include "pargs.h"
#else
#include "statutil.h"
#endif

/* Struct for a weighted grid and its associated information */
struct weighted_grid{
	real *weights; // Indexed [x][z][y] using pointer arithmetic (y is last index for y-axis maximization performance)
	int dimx, dimz, dimy; // Number of grid points in each dimension (this is number of grid cells + 1)
	real cell_width; // The length of one edge of a grid cell in trajectory space
	real minx, miny, minz; // The coordinates of the grid's origin in trajectory space
};

real tessellate_area(const char *traj_fname, const char *ndx_fname, int numcells, output_env_t *oenv);
/* Calculates the approximate surface area of a trajectory (or part of it specified by the given index file) by tessellating the coordinates in a 3D grid.
 * If ndx_fname is null, the whole trajectory will be included in the grid.
 * numcells is the number of grid cells to be created in the longest dimension.
 * Ideally, numcells should be lower than the number of trajectory points in the longest dimension.  
 * The approximated surface area is returned.
 */

void construct_grid(rvec **x, int nframes, int natoms, int numcells, struct weighted_grid *grid);
/* Memory is allocated for weights in grid and initialized to 0.
 * Call free_grid when done with grid.
 */

void load_grid(rvec **x, int nframes, int natoms, struct weighted_grid *grid);

void free_grid(struct weighted_grid *grid);

#endif // LL_TESSELLATION_H