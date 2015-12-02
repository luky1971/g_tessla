/*
 * Copyright Ahnaf Siddiqui and Sameer Varma.
 */

#ifndef LLT_TRI_H
#define LLT_TRI_H

#ifdef GRO_V5
#include "pargs.h"
#else
#include "statutil.h"
#endif

void llt_tri_area(const char *traj_fname, const char *ndx_fname, output_env_t *oenv);

#endif // LLT_TRI_H