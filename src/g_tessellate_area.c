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
        "g_tessellate_area calculates 3-d surface area using Delaunay tessellation. \n",
        "It reads in a trajectory file through the -f option (supported formats=xtc,trr,pdb). \n",
        "The set of points for tessellation, such as the coordinates of phosphorous atoms in a lipid bilayer, are specified using an index file by the -n option.\n",
        "Areas can be calculated individually for each frame in which case the output is dumped into an ASCII file specified by the -o option.\n\n",
        "This code can also be used for calculating the surface areas of lipid bilayers.\n", 
        "In such a calculation, the lipid bilayer normal is assumed to be parallel to the z-axis.\n",
        "This assumption is made to include in the surface area the space between the atoms lying at the periphery of the unit cell and the boundary of the unit cell.\n",
        "The correction is performed by inserting points at regular intervals along the edges of the simulation box.\n",
        "To use this correction, set the boolean -corr \n",
        "You can also set -espace X, ",
        "where X is the desired spacing in nanometers of the edge correction point intervals (default = 0.8).\n\n",
        "The -2d option will yield 2D projections on the XY plane - \n",
        "for a lipid bilayer perpendicular to the z-axis, the 2D projected area along with the -corr option ", 
        "will essentially yield the 2D area of the simulation cell.\n\n",
        "An alternative way to calculate lipid surface areas is to map the coordinates onto a weighted 3D grid \n", 
        "and tessellate the highest weight z-coordinates along the horizontal plane.\n",
        "The latter method is, however, still experimental and not supported. \n",
        "To use the experimental weighted grid method, set the -dense option.\n\n",       
        "The tessellated surface can be visualized using the -print option. The resulting .node and .ele files are numbered by frame \n",
        "and can be viewed by Jonathan R. Shewchuck's program showme\n",
        "(found here: https://www.cs.cmu.edu/~quake/showme.html)\n",
        "WARNING, the -print option produces a .node and .ele file for EVERY frame AND disables parallelization!\n",
        "(So don't be surprised when you come back hours later and see a hundred thousand new files in your current directory)\n\n",
        "If you build g_tessellate_area with OPENMP, you can set the number of threads to use with -nthreads X,\n",
        "where X is the number of threads to use. The default is to use the maximum number of cores available.\n"
    };

    const char *fnames[efT_NUMFILES];
    output_env_t oenv = NULL;

    int nthreads = -1;
    gmx_bool dense = FALSE;
    gmx_bool corr = FALSE;
    real espace = 0.8;
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
        {"-corr", FALSE, etBOOL, {&corr}, "correct triangulation area for periodic bounding"},
        {"-espace", FALSE, etREAL, {&espace}, "the spacing of the edge correction point intervals if using -corr (default = 0.8)"},
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

        unsigned long flags = ((int)corr * GTA_CORRECT) 
                            | ((int)a2D * GTA_2D) 
                            | ((int)print * GTA_PRINT);
        
        tessellate_area(fnames[efT_TRAJ], fnames[efT_NDX], &oenv, espace, nthreads, &areas, flags);

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
