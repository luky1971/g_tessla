/*
 * Copyright Ahnaf Siddiqui and Sameer Varma.
 */

#include "llt_tri.h"

#include "gkut_io.h"
#include "smalloc.h"
#include "triangle.h"

void llt_tri_area(const char *traj_fname, const char *ndx_fname, output_env_t *oenv) {
	rvec **pre_x, **x;
	int nframes, natoms;

	read_traj(traj_fname, &pre_x, &nframes, &natoms, oenv);

	// Filter trajectory by index file if present
	if(ndx_fname != NULL) {
		ndx_filter_traj(ndx_fname, &pre_x, &x, nframes, &natoms);
	}
	else {
		x = pre_x;
	}

	//

	// free memory
	for(int i = 0; i < nframes; ++i) {
		sfree(x[i]);
	}
	sfree(x);
}

//void f_llt_tri_area(rvec **x, )
