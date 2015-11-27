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

/* Struct for a tessellated grid and its associated information */
struct tessellated_grid{
	real *weights; // [dimx][dimy][dimz]
	int *heightmap; // [dimx][dimy]. Holds the z-index of the grid point with the maximum weight for each x-y column
	real *areas; // [dimx-1][dimy-1]. Holds the triangulated area of each grid cell
	int dimx, dimy, dimz; // Number of grid points in each dimension (this is number of grid cells + 1)
	real cell_width; // The length of one edge of a grid cell in trajectory space
	real minx, miny, minz; // The coordinates of the grid's origin in trajectory space
	int num_empty; // The number of grid cells that have zero-weight corner(s) (and are therefore excluded from tessellation)
	real surface_area; // The total triangulated/mesh surface area of the grid's heightmap
};

void llt_area(const char *traj_fname, const char *ndx_fname, real cell_width, output_env_t *oenv, struct tessellated_grid *grid);
/* Reads a trajectory file and then calculates approximate surface area (see the f_llt_area function below).
 * If ndx_fname is not null, only a selection within the trajectory will be included in the grid.
 */

void f_llt_area(rvec **x, int nframes, int natoms, real cell_width, struct tessellated_grid *grid);
/* Calculates the approximate surface area of a trajectory by tessellating the coordinates in a 3D grid.
 * Stores grid and area information in tessellated_grid *grid.
 * Call this if you have already read the trajectory (otherwise call llt_area above).
 * cell_width is the width of each grid cell. It should be high enough so that there's no empty grid cells inside the system of interest.
 * Memory is allocated for arrays in grid. Call free_grid when done.
 */

void construct_grid(rvec **x, int nframes, int natoms, real cell_width, struct tessellated_grid *grid);
/* Memory is allocated for arrays in grid and initialized to 0.
 * Call free_grid when done with grid.
 */

/* construct_grid must be called prior to calling the following functions */

void load_grid(rvec **x, int nframes, int natoms, struct tessellated_grid *grid);

void gen_heightmap(struct tessellated_grid *grid);

void tessellate_area(struct tessellated_grid *grid);
/* Tessellates the heightmap in the given grid and calculates the total area of the triangulated surface. 
 */ 

void print_grid(struct tessellated_grid *grid, const char *fname);
/* Prints the data in the given tessellated grid to the given file.
 */

void free_grid(struct tessellated_grid *grid);

#endif // LL_TESSELLATION_H