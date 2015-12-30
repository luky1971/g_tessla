/*
 * Copyright Ahnaf Siddiqui
 */

#include "gkut_io.h"

#ifdef GRO_V5
#include "index.h"
#include "trxio.h"
#endif
#include "smalloc.h"

void read_traj(const char *traj_fname, rvec ***x, matrix **box, int *nframes, int *natoms, output_env_t *oenv) {
	t_trxstatus *status = NULL;
	real t;
	int est_frames = FRAMESTEP;
	*nframes = 0;

	snew(*x, est_frames);
	snew(*box, est_frames);
	*natoms = read_first_x(*oenv, &status, traj_fname, &t, &((*x)[0]), (*box)[0]);

	do {
		++(*nframes);
		if(*nframes >= est_frames) {
			est_frames += FRAMESTEP;
			srenew(*x, est_frames);
			srenew(*box, est_frames);
		}
		snew((*x)[*nframes], *natoms);
	} while(read_next_x(*oenv, status, &t,
#ifndef GRO_V5 
		*natoms,
#endif
		(*x)[*nframes], (*box)[*nframes]));

	sfree((*x)[*nframes]); // Nothing was read to the last allocated frame
	close_trx(status);
}

void print_traj(rvec **x, int nframes, int natoms, const char *fname) {
	int fr, i;
	FILE *f = fopen(fname, "w");

	for(fr = 0; fr < nframes; ++fr) {
		fprintf(f, "\nFrame %d\n", fr);
		for(i = 0; i < natoms; ++i) {
			fprintf(f, "%d: %f %f %f\n", i, x[fr][i][XX], x[fr][i][YY], x[fr][i][ZZ]);
		}
	}

	fclose(f);
}

void ndx_get_indx(const char *ndx_fname, int numgroups, atom_id ***indx, int **isize) {
	char **grp_names;

	snew(*isize, numgroups);
	snew(*indx, numgroups);
	snew(grp_names, numgroups);

	rd_index(ndx_fname, numgroups, *isize, *indx, grp_names);
	sfree(grp_names);
}

void filter_vecs(atom_id *indx, int isize, rvec *pre_x, rvec **new_x) {
	snew(*new_x, isize);
	for(int i = 0; i < isize; ++i) {
		copy_rvec(pre_x[indx[i]], (*new_x)[i]);
	}
}

void ndx_filter_traj(const char *ndx_fname, rvec **pre_x, rvec ***new_x, int nframes, int *natoms) {
	const int NUMGROUPS = 1;
	int *isize;
	atom_id **indx;

	ndx_get_indx(ndx_fname, NUMGROUPS, &indx, &isize);

	*natoms = isize[0];
	sfree(isize);

	snew(*new_x, nframes);
	for(int i = 0; i < nframes; ++i) {
		filter_vecs(indx[0], *natoms, pre_x[i], &((*new_x)[i]));
	}

	sfree(indx[0]);
	sfree(indx);
}
