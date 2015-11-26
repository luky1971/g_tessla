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

#endif // GKUT_IO_H