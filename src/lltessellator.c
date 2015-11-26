/*
 * Copyright Ahnaf Siddiqui and Sameer Varma.
 */

#include "macros.h"

#include "gkut_log.h"

#include "ll_tessellation.h"

enum {efT_TRAJ, efT_NDX, efT_NUMFILES};

int main(int argc, char *argv[]) {
	const char *desc[] = {
		"lltessellator reads a trajectory file and tessellates its coordinates in a fixed 3D grid,",
		" calculating the surface area of the tessellation.\n"
	};
	const char *fnames[efT_NUMFILES];
	output_env_t oenv = NULL;
	real cell_width = 0.1;

	init_log("llt.log", argv[0]);

	t_filenm fnm[] = {
		{efTRX, "-f", "traj.xtc", ffREAD},
		{efNDX, "-n", "index.ndx", ffOPTRD}
	};

	t_pargs pa[] = {
		{"-width", FALSE, etREAL, {&cell_width}, "width of each grid cell."}
	};

	parse_common_args(&argc, argv, 0, efT_NUMFILES, fnm, asize(pa), pa, asize(desc), desc, 0, NULL, &oenv);

	fnames[efT_TRAJ] = opt2fn("-f", efT_NUMFILES, fnm);
	fnames[efT_NDX] = opt2fn_null("-n", efT_NUMFILES, fnm);

	tessellate_area(fnames[efT_TRAJ], fnames[efT_NDX], cell_width, &oenv);

	print_log("Done.\n");
	close_log();

	return 0;
}
