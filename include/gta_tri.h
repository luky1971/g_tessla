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

#ifndef GTA_TRI_H
#define GTA_TRI_H

#include "vec.h"
#ifdef GRO_V5
#include "pargs.h"
#else
#include "statutil.h"
#endif


// Flags
enum {
    GTA_CORRECT = 1, // Correct areas for periodic bounding conditions
    GTA_2D = 2, // Calculate 2D surface area as well
    GTA_PRINT = 4, // Print triangle data that can be visualized using, for example, the 'showme' program
};

// Struct for area output data.
// These are total surface area, divide a given area by natoms to get area per particle.
struct tri_area {
    real *area; // Triangulated 3D areas indexed by [frame #]. *area are corrected areas for periodic bounds if GTA_CORRECT was used.
    real *area2D; // Triangulated 2D areas indexed by [frame #]. NULL if GTA_2D not set.
    real *area2Dbox; // 2D areas of box for each frame.
    int natoms, nframes; // Number of atoms and number of frames, respectively, that were triangulated.
};


void tessellate_area(const char *traj_fname, 
                     const char *ndx_fname, 
                     output_env_t *oenv, 
                     real espace, 
                     int nthreads, 
                     struct tri_area *areas, 
                     unsigned char flags);
/* Reads a trajectory file and tessellates all of its frames.
 * If ndx_fname is not null, only a selection within the trajectory will be tessellated.
 * output_env_t *oenv is needed for reading trajectory files.
 * You can initialize one using output_env_init() in Gromacs's oenv.h.
 * Calls the delaunay_tessellate function below.
 */

void delaunay_tessellate(rvec **x, 
                         matrix *box, 
                         real espace, 
                         int nthreads, 
                         struct tri_area *areas, 
                         unsigned char flags);
/* Tesssellates all of the frames in the given trajectory using delaunay triangulation.
 * espace is the spacing of the edge correction point intervals if using the GTA_CORRECT flag.
 * nthreads is the number of threads to be used if ensemble_comp was built using openmp.
 * nthreads <= 0 will use all available threads.
 * Memory is allocated for arrays in the tri_area struct. Call free_tri_area when done.
 * See above for flags.
 */

void delaunay_surface_area(const rvec *x, 
                           matrix box, 
                           int natoms, 
                           unsigned char flags, 
                           real *a2D, 
                           real *a3D);
/* Tessellates the given array of coordinates using delaunay triangulation 
 * and calculates 2D and 3D area, stored in a2D and a3D.
 * a2D and/or a3D can be NULL.
 * See above for flags.
 */

void print_areas(const char *fname, const struct tri_area *areas);
/* Formats and prints the data in a tri_area struct to an output file.
 */

void free_tri_area(struct tri_area *areas);
/* Frees the dynamic memory in a tri_area struct.
 */


/* Calculates the area of the triangle formed by three points in 3D space.
 */
static inline real area_tri(const rvec a, 
                            const rvec b, 
                            const rvec c) {
    rvec ab, ac, cpr;
    rvec_sub(b, a, ab);
    rvec_sub(c, a, ac);
    cprod(ab, ac, cpr);
    return norm(cpr) / 2.0;
}


#endif // GTA_TRI_H