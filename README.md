### lltessellator

lltessellator (leaflet tessellator) reads a trajectory file, tessellates its coordinates,
and calculates 3D surface area. 

The trajectory should represent a planar fluctuating surface such as a lipid leaflet.
It is specified by the -f option (supported formats=xtc,trr,pdb), and can be filtered using an index file
by the -n option. An output filename can be specified with the -o option. 

lltessellator can either calculate every frame's 3D surface area using delaunay triangulation on the horizontal plane,
or it can load the points from all frames into a weighted 3D grid and tessellate the grid
based on highest weight z-coordinates along the horizontal plane.
This latter method is experimental and not supported; delaunay triangulation, the default, is more accurate.
To use the experimental weighted grid method, set the -dense option. 

When using delaunay triangulation, the surface area can be corrected for periodic bounding conditions.
The correction is performed by inserting points at a given interval along the edges of the simulation box
with z-coordinates that are the weighted average of the two points closest to the two edges in that interval.
Points are inserted at the corners of the box with the same z coordinate, which is the average of the four points
closest to the four corners of the box.
To use this correction method for periodic bounding conditions, set the option -corr X
where X is the desired spacing of the edge correction point intervals. 

The -2 option for calculating 2D area from triangulation is mostly for testing purposes,
as this should be equivalent to the 2D box area which is calculated and output anyways. 

If you set the -print option for delaunay triangulation, the resulting .node and .ele files are numbered by frame 
and can be viewed by Jonathan R. Shewchuck's program showme (found here: https://www.cs.cmu.edu/~quake/showme.html).
BE WARNED, the -print option produces a .node and .ele file for EVERY frame AND disables parallelization!
(So don't be surprised when you come back hours later and see a hundred thousand new files in your current directory)

By default, lltessellator is parallelized with OpenMP (see installation instructions below).
To disable parallelization at runtime, add the option -nopar.

### INSTALLATION

The following instructions are for unix-based operating systems such as OSX and Linux.
Windows is not currently supported.

1. Install Gromacs version 4.5.x or later from http://www.gromacs.org.

2. `git clone` or otherwise obtain and `cd` to the 'lltessellator' repository.

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
you may have to add lltessellator to your PATH. 

If you are not using gcc, you will also need to set `CC` and `CXX`
to your C compiler and C++ compiler commands respectively.

If you want to build without OpenMP, set `PARALLEL=0`. You can also add compilation flags by setting `CFLAGS`.

### Copyright 
(c) 2016 Ahnaf Siddiqui and Sameer Varma 

lltessellator's implementation of Delaunay triangulation is based on the algorithms presented by 

Lee, D.T. and Schachter, B.J. Two algorithms for constructing a Delaunay triangulation.
International Journal of Computer & Information Sciences 1980;9(3):219-242. 

and 

Guibas, L. and Stolfi, J. Primitives for the manipulation of general subdivisions and the computation of Voronoi.
ACM Trans. Graph. 1985;4(2):74-123. 

lltessellator uses exact arithmetic routines and geometric predicates provided by 

Shewchuk, J.R. 1996. Routines for Arbitrary Precision Floating-point Arithmetic and Fast Robust Geometric Predicates.

lltessellator uses the GROMACS molecular simulation package API.  
Copyright (c) 1991-2000, University of Groningen, The Netherlands.  
Copyright (c) 2001-2004, The GROMACS development team.  
Copyright (c) 2013,2014, by the GROMACS development team,  
led by Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,  
and including many others, as listed at http://www.gromacs.org.