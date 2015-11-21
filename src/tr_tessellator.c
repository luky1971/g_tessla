/*
 * This program uses the GROMACS molecular simulation package API.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.

 * tr_tessellator reads a trajectory file and tessellates its coordinates in a fixed 3D grid, calculating the surface area of the tessellation.
 * Written by Ahnaf Siddiqui and Sameer Varma.
 */

#include "macros.h"
#ifdef GRO_V5
#include "pargs.h"
#else
#include "statutil.h"
#endif

enum {efT_TRAJ, efT_NDX, efT_NUMFILES};

int main(int argc, char *argv[]) {
	const char *desc[] = {
		"tr_tessellator reads a trajectory file and tessellates its coordinates in a fixed 3D grid,",
		" calculating the surface area of the tessellation.\n"
	};
	const char *fnames[efT_NUMFILES];
	output_env_t oenv = NULL;
	real width = 0.1; // width of each cell in the grid

	t_filenm fnm[] = {
		{efTRX, "-f", "traj.xtc", ffREAD},
		{efNDX, "-n", "index.ndx", ffOPTRD}
	};

	t_pargs pa[] = {
		{"-w", FALSE, etREAL, {&width}, "Width of each cell in the tessellation grid"}
	};

	parse_common_args(&argc, argv, 0, efT_NUMFILES, fnm, asize(pa), pa, asize(desc), desc, 0, NULL, &oenv);

	fnames[efT_TRAJ] = opt2fn("-f", efT_NUMFILES, fnm);
	fnames[efT_NDX] = opt2fn_null("-n", efT_NUMFILES, fnm);

	return 0;
}
