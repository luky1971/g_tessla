/*
 * Copyright Ahnaf Siddiqui
 */

#ifndef GKUT_IO_H
#define GKUT_IO_H

#include "vec.h"
#ifdef GRO_V5
#include "pargs.h"
#else
#include "statutil.h"
#endif

#define FRAMESTEP 500 // The number of new frames by which to reallocate an array of length # trajectory frames

void read_traj(const char *traj_fname, rvec ***x, int *nframes, int *natoms, output_env_t *oenv);
/* Reads a trajectory file. rvec **x is position coordinates indexed x[frame #][atom #]
 * 2D memory is allocated for x.
 */

void ndx_filter_traj(const char *ndx_fname, rvec ***pre_x, rvec ***new_x, int nframes, int *natoms);
/* Creates a new trajectory with only the coordinates from the old trajectory (pre_x) that are specified by the index file.
 * 2D memory is freed for pre_x and allocated for x.
 */

#endif // GKUT_IO_H