/*
 * Copyright Ahnaf Siddiqui and Sameer Varma.
 */

#ifndef TRT_TESSELLATION_H
#define TRT_TESSELLATION_H

real tr_tessellate_area(const char *traj_fname, const char *ndx_fname);
/* Calculates the approximate surface area of a trajectory (or part of it specified by the given index file) by tessellating the coordinates in a 3D grid.
 * If ndx_fname is null, the whole trajectory will be included in the grid.
 * The approximated surface area is returned.
 */

#endif // TRT_TESSELLATION_H