/*
 * Copyright Ahnaf Siddiqui and Sameer Varma.
 */

#ifndef LLT_TRI_H
#define LLT_TRI_H

#include "triangle.h"
#include "vec.h"
#ifdef GRO_V5
#include "pargs.h"
#else
#include "statutil.h"
#endif

void llt_tri_area(const char *traj_fname, const char *ndx_fname, output_env_t *oenv, 
	real **areas, int *nframes, int *natoms);

real tri_surface_area(rvec *x, int natoms);

void print_areas(const char *fname, real *areas, int nframes, int natoms);

void print_trifiles(struct triangulateio *tio, const char *node_name, const char *ele_name);


inline real area_tri(rvec a, rvec b, rvec c) {
	rvec ab, ac, cpr;
	rvec_sub(b, a, ab);
	rvec_sub(c, a, ac);
	cprod(ab, ac, cpr);
	return norm(cpr) / 2.0;
}


#endif // LLT_TRI_H