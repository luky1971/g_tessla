/*
 * Copyright Ahnaf Siddiqui
 */

#include "gkut_io.h"

#include "smalloc.h"

void read_traj(const char *traj_fname, rvec ***x, int *nframes, int *natoms, output_env_t *oenv) {
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

void ndx_filter_traj(const char *ndx_fname, rvec ***pre_x, rvec ***new_x, int nframes, int *natoms) {
	const int NUMGROUPS = 1;
	int *isize;
	atom_id **indx;
	char **grp_names;

	snew(isize, NUMGROUPS);
	snew(indx, NUMGROUPS);
	snew(grp_names, NUMGROUPS);

	rd_index(ndx_fname, NUMGROUPS, isize, indx, grp_names);
	sfree(grp_names);

	*natoms = isize[0];
	sfree(isize);

	snew(*new_x, nframes);
	for(int i = 0; i < nframes; ++i) {
		snew((*new_x)[i], *natoms);
		for(int j = 0; j < *natoms; ++j) {
			copy_rvec((*pre_x)[i][indx[0][j]], (*new_x)[i][j]);
		}
	}

	// free pre (unfiltered) trajectory's memory
	sfree(indx[0]);
	sfree(indx);

	for(int i = 0; i < nframes; ++i) {
		sfree((*pre_x)[i]);
	}
	sfree(*pre_x);
}
