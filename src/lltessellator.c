/*
 * Copyright Ahnaf Siddiqui and Sameer Varma.
 */

#ifdef LLT_BENCH
#include <time.h>
#endif
#include "macros.h"
#include "smalloc.h"

#include "gkut_log.h"

#include "llt_grid.h"
#include "llt_tri.h"

enum {efT_TRAJ, efT_NDX, efT_OUTDAT, efT_NUMFILES};

int main(int argc, char *argv[]) {
#ifdef LLT_BENCH
	clock_t start = clock();
#endif

	const char *desc[] = {
		"lltessellator reads a trajectory file and tessellates its coordinates.\n",
		"It can either tessellate the points in each frame using delaunay triangulation and calculate average surface area,\n",
		"or it can load the points from all frames into a weighted grid and tessellate the grid based on density.\n",
		"Set the -dense option to use the latter method.\n"
	};
	const char *fnames[efT_NUMFILES];
	output_env_t oenv = NULL;

	gmx_bool dense = FALSE;
	gmx_bool correct = FALSE;
	gmx_bool a2D = FALSE;
	gmx_bool print = FALSE;
	real cell_width = 0.1;
	gmx_bool linear = FALSE;

	init_log("llt.log", argc, argv);

	t_filenm fnm[] = {
		{efTRX, "-f", "traj.xtc", ffREAD},
		{efNDX, "-n", "index.ndx", ffOPTRD},
		{efDAT, "-o", "tessellation.dat", ffWRITE}
	};

	t_pargs pa[] = {
		{"-dense", FALSE, etBOOL, {&dense}, "use weighted-grid tessellation instead of frame-by-frame triangulation"},
		{"-correct", FALSE, etBOOL, {&correct}, "correct triangulation area for periodic bounding conditions (for delaunay triangulation)"},
		{"-2", FALSE, etBOOL, {&a2D}, "calculate 2D surface area from triangulation"},
		{"-print", FALSE, etBOOL, {&print}, "save triangles to .node and .ele files (for delaunay triangulation)"},
		{"-width", FALSE, etREAL, {&cell_width}, "width of each grid cell if using -dense"},
		{"-lin", FALSE, etBOOL, {&linear}, "use distance instead of distance squared for weighing if using -dense"}
	};

	parse_common_args(&argc, argv, 0, efT_NUMFILES, fnm, asize(pa), pa, asize(desc), desc, 0, NULL, &oenv);

	fnames[efT_TRAJ] = opt2fn("-f", efT_NUMFILES, fnm);
	fnames[efT_NDX] = opt2fn_null("-n", efT_NUMFILES, fnm);
	fnames[efT_OUTDAT] = opt2fn("-o", efT_NUMFILES, fnm);

	if(dense) {
		real (*fweight)(rvec, rvec) = linear ? weight_dist : weight_dist2;
		struct tessellated_grid grid;

		llt_grid_area(fnames[efT_TRAJ], fnames[efT_NDX], cell_width, fweight, &oenv, &grid);

		if(grid.num_empty > 0) {
			print_log("\n\nWARNING: %d grid cell(s) have empty corner(s).\n"
				"If there are gaps in your system's area of the grid, increase grid spacing with the -width option.\n"
				"The current width is %f.\n", grid.num_empty, grid.cell_width);
		}

		print_log("Tessellated surface area per particle: %f\n", grid.area_per_particle);

		print_grid(&grid, fnames[efT_OUTDAT]);

		free_grid(&grid);
	}
	else {
		struct tri_area areas;

		unsigned long flags = ((int)correct * LLT_CORRECT) | ((int)a2D * LLT_2D) | ((int)print * LLT_PRINT);
		
		llt_tri_area(fnames[efT_TRAJ], fnames[efT_NDX], &oenv, &areas, flags);

		print_areas(fnames[efT_OUTDAT], &areas);

		free_tri_area(&areas);
	}

#ifdef LLT_BENCH
	clock_t clocks = clock() - start;
	print_log("lltessellator took %d clocks, %f seconds.\n", clocks, (float)clocks/CLOCKS_PER_SEC);
#endif

	close_log();

	return 0;
}
