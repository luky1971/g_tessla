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

#ifndef LLT_TRI_H
#define LLT_TRI_H

#include "vec.h"
#ifdef GRO_V5
#include "pargs.h"
#else
#include "statutil.h"
#endif


// Flags
enum {
    LLT_NOPAR = 1, // Prevents multithreading
    LLT_CORRECT = 2, // Correct areas for periodic bounding conditions
    LLT_2D = 4, // Calculate 2D surface area as well
    LLT_PRINT = 8 // Print triangle data that can be visualized using, for example, the 'showme' program
};

// Struct for area output data.
// These are total surface area, divide a given area by natoms to get area per particle.
struct tri_area {
    real *area; // Triangulated 3D areas indexed by [frame #]. *area are corrected areas for periodic bounds if LLT_CORRECT was used.
    real *area2D; // Triangulated 2D areas indexed by [frame #]. NULL if LLT_2D not set.
    real *area2Dbox; // 2D areas of box for each frame.
    int natoms, nframes; // Number of atoms and number of frames, respectively, that were triangulated.
};


void llt_delaunay_area( const char *traj_fname, 
                        const char *ndx_fname, 
                        output_env_t *oenv, 
                        real corr, 
                        struct tri_area *areas, 
                        unsigned char flags);
/* Reads a trajectory file and tessellates all of its frames using delaunay triangulation.
 * If ndx_fname is not null, only a selection within the trajectory will be tessellated.
 * output_env_t *oenv is needed for reading trajectory files.
 * You can initialize one using output_env_init() in Gromacs's oenv.h.
 * Memory is allocated for arrays in the tri_area struct. Call free_tri_area when done.
 * See above for flags.
 */

void delaunay_surface_area( const rvec *x, 
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


#endif // LLT_TRI_H