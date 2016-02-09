/*
 * Copyright Ahnaf Siddiqui and Sameer Varma.
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


void llt_delaunay_area(	const char *traj_fname, 
						const char *ndx_fname, 
						output_env_t *oenv, 
						real corr, 
						struct tri_area *areas, 
						unsigned char flags);

void delaunay_surface_area(	const rvec *x, 
							int natoms, 
							unsigned char flags,
							real *a2D,
							real *a3D);
/* a2D and/or a3D can be NULL.
 */

void print_areas(const char *fname, const struct tri_area *areas);

void free_tri_area(struct tri_area *areas);


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