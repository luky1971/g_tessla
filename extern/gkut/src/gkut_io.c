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
