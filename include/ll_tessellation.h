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
	real *weights; // Indexed [x][y][z] using pointer arithmetic
	int dimx, dimy, dimz; // Number of grid points in each dimension (this is number of grid cells + 1)
	real cell_width; // The length of one edge of a grid cell in trajectory space
	real minx, miny, minz; // The coordinates of the grid's origin in trajectory space
};

real tessellate_area(const char *traj_fname, const char *ndx_fname, real cell_width, output_env_t *oenv);
/* Reads a trajectory file and then calculates approximate surface area (see the tessellate_area function below).
 * If ndx_fname is null, the whole trajectory will be included in the grid.
 */

real f_tessellate_area(rvec **x, int nframes, int natoms, real cell_width);
/* Calculates the approximate surface area of a trajectory by tessellating the coordinates in a 3D grid.
 * Call this if you have already read the trajectory.
 * cell_width is the width of each grid cell. It should be high enough so that there's no empty grid cells.
 * The approximated surface area is returned.
 */

void construct_grid(rvec **x, int nframes, int natoms, real cell_width, struct weighted_grid *grid);
/* Memory is allocated for weights in grid and initialized to 0.
 * Call free_grid when done with grid.
 */

void load_grid(rvec **x, int nframes, int natoms, struct weighted_grid *grid);

void free_grid(struct weighted_grid *grid);

#endif // LL_TESSELLATION_H