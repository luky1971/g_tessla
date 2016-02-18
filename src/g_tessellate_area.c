/*
 * Copyright 2016 Ahnaf Siddiqui and Sameer Varma
 *
 * This program uses the GROMACS molecular simulation package API.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed at http://www.gromacs.org.
 */

#ifdef GTA_BENCH
#include <time.h>
#endif
#include "macros.h"
#include "smalloc.h"

#include "gkut_log.h"

#include "gta_grid.h"
#include "gta_tri.h"

#define CORR_EPS 1e-12

enum {efT_TRAJ, efT_NDX, efT_OUTDAT, efT_NUMFILES};

int main(int argc, char *argv[]) {
#ifdef GTA_BENCH
    clock_t start = clock();
#endif

    const char *desc[] = {
        "g_tessellate_area reads a trajectory file, tessellates its coordinates, \n",
        "and calculates 3D surface area.\n\n",
        "The trajectory should represent a planar fluctuating surface such as a lipid leaflet.\n",
        "It is specified by the -f option (supported formats=xtc,trr,pdb), and can be filtered using an index file \n",
        "by the -n option. An output filename can be specified with the -o option.\n\n",
        "g_tessellate_area can either calculate every frame's 3D surface area using delaunay triangulation on the horizontal plane, \n",
        "or it can load the points from all frames into a weighted 3D grid and tessellate the grid \n",
        "based on highest weight z-coordinates along the horizontal plane.\n",
        "This latter method is experimental and not supported; delaunay triangulation, the default, is more accurate.\n",
        "To use the experimental weighted grid method, set the -dense option.\n\n",
        "When using delaunay triangulation, the surface area can be corrected for periodic bounding conditions.\n",
        "The correction is performed by inserting points at a given interval along the edges of the simulation box \n",
        "with z-coordinates that are the weighted average of the two points closest to the two edges in that interval.\n",
        "Points are inserted at the corners of the box with the same z coordinate, which is the average of the four points \n",
        "closest to the four corners of the box.\n",
        "To use this correction method for periodic bounding conditions, set the option -corr X \n",
        "where X is the desired spacing of the edge correction point intervals.\n\n",
        "The -2d option for calculating 2D area from triangulation is mostly for testing purposes, \n",
        "as this should be equivalent to the 2D box area which is calculated and output anyways.\n\n",
        "If you set the -print option for delaunay triangulation, the resulting .node and .ele files are numbered by frame \n",
        "and can be viewed by Jonathan R. Shewchuck's program showme\n",
        "(found here: https://www.cs.cmu.edu/~quake/showme.html)\n",
        "BE WARNED, the -print option produces a .node and .ele file for EVERY frame AND disables parallelization!\n",
        "(So don't be surprised when you come back hours later and see a hundred thousand new files in your current directory)\n\n",
        "If g_tessellate_area was built with OPENMP, you can set the number of threads to use with -nthreads X,\n",
        "where X is the number of threads to use. The default behavior is to use the maximum number of cores available.\n"
    };
    const char *fnames[efT_NUMFILES];
    output_env_t oenv = NULL;

    int nthreads = -1;
    gmx_bool dense = FALSE;
    real corr = 0.0;
    gmx_bool a2D = FALSE;
    gmx_bool print = FALSE;
    real cell_width = 0.1;
    gmx_bool linear = FALSE;

    init_log("gta.log", argc, argv);

    t_filenm fnm[] = {
        {efTRX, "-f", "traj.xtc", ffREAD},
        {efNDX, "-n", "index.ndx", ffOPTRD},
        {efDAT, "-o", "tessellated_areas.dat", ffWRITE}
    };

    t_pargs pa[] = {
        {"-nthreads", FALSE, etINT, {&nthreads}, "set the number of parallel threads to use (default is max available)"}, 
        {"-dense", FALSE, etBOOL, {&dense}, "use weighted-grid tessellation instead of frame-by-frame delaunay triangulation"},
        {"-corr", FALSE, etREAL, {&corr}, "correct triangulation area for periodic bounding by given edge spacing (see readme)"},
        {"-2d", FALSE, etBOOL, {&a2D}, "calculate 2D surface area from delaunay triangulation"},
        {"-print", FALSE, etBOOL, {&print}, "BE CAREFUL (see readme); save delaunay triangles to .node and .ele files"},
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

        gta_grid_area(fnames[efT_TRAJ], fnames[efT_NDX], cell_width, fweight, &oenv, &grid);

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
        if(print)   nthreads = 1;

        struct tri_area areas;

        unsigned long flags = ((int)(corr > CORR_EPS) * GTA_CORRECT) 
                            | ((int)a2D * GTA_2D) 
                            | ((int)print * GTA_PRINT);
        
        gta_delaunay_area(fnames[efT_TRAJ], fnames[efT_NDX], &oenv, corr, nthreads, &areas, flags);

        print_areas(fnames[efT_OUTDAT], &areas);

        free_tri_area(&areas);
    }

#ifdef GTA_BENCH
    clock_t clocks = clock() - start;
    print_log("g_tessellate_area took %d clocks, %f seconds.\n", clocks, (float)clocks/CLOCKS_PER_SEC);
#endif

    close_log();

    return 0;
}
