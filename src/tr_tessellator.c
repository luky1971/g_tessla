/*
 * Copyright Ahnaf Siddiqui and Sameer Varma.
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
