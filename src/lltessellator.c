/*
 * Copyright Ahnaf Siddiqui and Sameer Varma.
 */

#include "macros.h"

#include "gkut_log.h"

#include "ll_tessellation.h"

enum {efT_TRAJ, efT_NDX, efT_OUTDAT, efT_NUMFILES};

int main(int argc, char *argv[]) {
	const char *desc[] = {
		"lltessellator reads a trajectory file and tessellates its coordinates in a fixed 3D grid,",
		" calculating the surface area of the tessellation.\n"
	};
	const char *fnames[efT_NUMFILES];
	output_env_t oenv = NULL;
	real cell_width = 0.1;
	real (*fweight)(rvec, rvec) = weight_dist2;
	gmx_bool linear = FALSE;

	struct tessellated_grid grid;

	init_log("llt.log", argv[0]);

	t_filenm fnm[] = {
		{efTRX, "-f", "traj.xtc", ffREAD},
		{efNDX, "-n", "index.ndx", ffOPTRD},
		{efDAT, "-o", "grid.dat", ffWRITE}
	};

	t_pargs pa[] = {
		{"-width", FALSE, etREAL, {&cell_width}, "width of each grid cell"},
		{"-lin", FALSE, etBOOL, {&linear}, "set this to use distance instead of distance squared for weighing"}
	};

	parse_common_args(&argc, argv, 0, efT_NUMFILES, fnm, asize(pa), pa, asize(desc), desc, 0, NULL, &oenv);

	fnames[efT_TRAJ] = opt2fn("-f", efT_NUMFILES, fnm);
	fnames[efT_NDX] = opt2fn_null("-n", efT_NUMFILES, fnm);
	fnames[efT_OUTDAT] = opt2fn("-o", efT_NUMFILES, fnm);

	if(linear)	fweight = weight_dist;

	llt_area(fnames[efT_TRAJ], fnames[efT_NDX], cell_width, fweight, &oenv, &grid);

	if(grid.num_empty > 0) {
		print_log("\n\nWARNING: %d grid cell(s) have empty corner(s).\n"
			"If there are gaps in your system's area of the grid, increase grid spacing with the -width option.\n"
			"The current width is %f.\n", grid.num_empty, grid.cell_width);
	}

	print_grid(&grid, fnames[efT_OUTDAT]);

	print_log("Tessellation data saved to %s\n", fnames[efT_OUTDAT]);
	print_log("Tessellated surface area: %f\n", grid.surface_area);

	close_log();

	return 0;
}
