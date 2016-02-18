/*
 * Copyright 2016 Ahnaf Siddiqui and Sameer Varma
 *
 * This program uses the GROMACS molecular simulation package API.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed at http://www.gromacs.org.
 */

#ifndef GTA_GRID_H
#define GTA_GRID_H

#include "vec.h"
#ifdef GRO_V5
#include "pargs.h"
#else
#include "statutil.h"
#endif

/* Struct for a weighted 3D grid and its associated information.
 * Note: some of the arrays in this struct are multidimensional as indicated by the comments,
 * but are stored as single pointers and addressed using pointer arithmetic.
 */
struct tessellated_grid {
    real *weights; // [dimx][dimy][dimz]
    int *heightmap; // [dimx][dimy]. Holds the z-index of the grid point with the maximum weight for each x-y column
    real *areas; // [dimx-1][dimy-1]. Holds the triangulated area of each grid cell
    int dimx, dimy, dimz; // Number of grid points in each dimension (this is number of grid cells + 1)
    real cell_width; // The length of one edge of a grid cell in trajectory space
    real minx, miny, minz; // The coordinates of the grid's origin in trajectory space
    int num_empty; // The number of grid cells that have zero-weight corner(s) (and are therefore excluded from tessellation)
    real surface_area; // The total triangulated/mesh surface area of the grid's heightmap
    real area_per_particle; // The mesh surface area per trajectory particle
};


/* Weight functions */

real weight_dist(rvec traj_point, rvec grid_point);
/* Assigns weight based on distance. Closer distance = higher weight.
 */

real weight_dist2(rvec traj_point, rvec grid_point);
/* Assigns weight based on distance squared. Closer distance squared = higher weight.
 */

/***/


void gta_grid_area(const char *traj_fname, const char *ndx_fname, 
    real cell_width, real (*fweight)(rvec, rvec), output_env_t *oenv, struct tessellated_grid *grid);
/* Reads a trajectory file and then calculates approximate surface area (see the f_gta_grid_area function below).
 * If ndx_fname is not null, only a selection within the trajectory will be included in the grid.
 * output_env_t *oenv is needed for reading trajectory files.
 * You can initialize one using output_env_init() in Gromacs's oenv.h.
 * Memory is allocated for arrays in grid. Call free_grid when done.
 */

void f_gta_grid_area(rvec **x, int nframes, int natoms, 
    real cell_width, real (*fweight)(rvec, rvec), struct tessellated_grid *grid);
/* Calculates the approximate surface area of a trajectory by tessellating the coordinates in a 3D grid.
 * Stores grid and area information in tessellated_grid *grid.
 * Call this if you have already read the trajectory (otherwise call gta_grid_area above).
 * cell_width is the width of each grid cell. It should be high enough so that there's no gaps (empty cells) within the system of interest.
 * fweight is the function that will be used for calculating the weight of each grid point - trajectory point pair.
 * You can use one of the weight functions above for fweight or your own weight function.
 * Memory is allocated for arrays in grid. Call free_grid when done.
 */

void construct_grid(rvec **x, int nframes, int natoms, real cell_width, struct tessellated_grid *grid);
/* Memory is allocated for arrays in grid and initialized to 0.
 * Call free_grid when done with grid.
 */

/* construct_grid must be called prior to calling the following functions */

void load_grid(rvec **x, int nframes, int natoms, real (*fweight)(rvec, rvec), struct tessellated_grid *grid);
/* Loads the given grid with weights based on the given trajectory.
 * Uses fweight to calculate the weight of each grid point - trajectory point pair.
 * You can use one of the weight functions above for fweight.
 */

void gen_heightmap(struct tessellated_grid *grid);
/* Finds the z-index of the grid point with the maximum weight for each x-y column in the grid.
 * This data is stored in grid->heightmap
 */

void tessellate_grid(struct tessellated_grid *grid);
/* Tessellates the heightmap in the given grid and calculates the total area of the triangulated surface. 
 */

void print_grid(struct tessellated_grid *grid, const char *fname);
/* Prints the data in the given tessellated grid to the given file.
 */

void free_grid(struct tessellated_grid *grid);

#endif // GTA_GRID_H