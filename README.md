### g_tessellate_area

g_tessellate_area calculates 3-d surface area using Delaunay tessellation.

It reads in a trajectory file through the -f option (supported formats=xtc,trr,pdb).
The set of points for tessellation, such as the coordinates of phosphorous atoms in a lipid bilayer, are specified using an index file by the -n option.
Areas can be calculated individually for each frame in which case the output is dumped into an ASCII file specified by the -o option. 

This code can also be used for calculating the surface areas of lipid bilayers.
In such a calculation, the lipid bilayer normal is assumed to be parallel to the z-axis.
This assumption is made to include in the surface area the space between the atoms lying at the periphery of the unit cell and the boundary of the unit cell.
The correction is performed by inserting points at regular intervals along the edges of the simulation box.
To use this correction, set the boolean -corr.
You can also set -espace X, where X is the desired spacing in nanometers of the edge correction point intervals (default = 0.8).

The -2d option will yield 2D projections on the XY plane - for a lipid bilayer perpendicular to the z-axis, the 2D projected area along with the -corr option will essentially yield the 2D area of the simulation cell.

An alternative way to calculate lipid surface areas is to map the coordinates onto a weighted 3D grid, and tessellate the highest weight z-coordinates along the horizontal plane. The latter method is, however, still experimental and not supported. To use the experimental weighted grid method, set the -dense option.

The tessellated surface can be visualized using the -print option. The resulting .node and .ele files are numbered by frame and can be viewed by Jonathan R. Shewchuck's program showme (found here: https://www.cs.cmu.edu/~quake/showme.html).
WARNING, the -print option produces a .node and .ele file for EVERY frame AND disables parallelization!
(So don't be surprised when you come back hours later and see a hundred thousand new files in your current directory)

If you build g_tessellate_area with OPENMP, you can set the number of threads to use with -nthreads X, where X is the number of threads to use. The default is to use the maximum number of cores available.

### INSTALLATION

The following instructions are for unix-based operating systems such as OSX and Linux.
Windows is not currently supported.

1. Install Gromacs version 4.5.x or later from http://www.gromacs.org.

2. `git clone` or otherwise obtain and `cd` to the 'g_tessellate_area' repository.

3. Run `sudo make install` with the necessary arguments for your environment (see below).

If you do not have Gromacs version 5.x installed, you will need to set the makefile's `VGRO`
variable to the root Gromacs version number.
If Gromacs is installed in a non-default directory (ie not in /usr/local/gromacs)
then you will have to set the `GROMACS` variable to the Gromacs installation directory
that contains the 'include' and 'lib' folders.
For example, if you are running Gromacs 4.5.3 installed in /home/user/tools/gromacs-4.5.3,
then you would run the following command:

`sudo make install VGRO=4 GROMACS=/home/user/tools/gromacs-4.5.3`

If you must run `make install` without sudo privileges, you will need to set the `INSTALL`
variable to a path that you can write to. 
The default install path is /usr/local/bin. Depending on your system and chosen installation directory,
you may have to add g_tessellate_area to your PATH. 

If you are not using gcc, you will also need to set `CC` and `CXX`
to your C compiler and C++ compiler commands respectively.

If you want to build without OpenMP, set `PARALLEL=0`. You can also add compilation flags by setting `CFLAGS`.

### Copyright 
(c) 2016 Ahnaf Siddiqui and Sameer Varma 

g_tessellate_area's implementation of Delaunay triangulation is based on the algorithms presented by 

Lee, D.T. and Schachter, B.J. Two algorithms for constructing a Delaunay triangulation.
International Journal of Computer & Information Sciences 1980;9(3):219-242. 

and 

Guibas, L. and Stolfi, J. Primitives for the manipulation of general subdivisions and the computation of Voronoi.
ACM Trans. Graph. 1985;4(2):74-123. 

g_tessellate_area uses exact arithmetic routines and geometric predicates provided by 

Shewchuk, J.R. 1996. Routines for Arbitrary Precision Floating-point Arithmetic and Fast Robust Geometric Predicates.

g_tessellate_area uses the GROMACS molecular simulation package API.  
Copyright (c) 1991-2000, University of Groningen, The Netherlands.  
Copyright (c) 2001-2004, The GROMACS development team.  
Copyright (c) 2013,2014, by the GROMACS development team,  
led by Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,  
and including many others, as listed at http://www.gromacs.org.