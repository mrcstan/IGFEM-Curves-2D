IGFEM-Curves-2D ========

This project calculates the temperature distribution in a 2D rectangular domain
given a channel network embedded in the domain. With the file describing the
channel in directory ChannelFiles (examples are files with ".channel" extension)
, you run "main.m" from the current directory

Features
--------

- Self-contained: does not need Abaqus to generate input file
- Efficient (typically a couple of minutes for each analysis)

Installation
------------

Two sets of C codes need to be compiled and interfaced with MATLAB. The
interfacing is done via Mex functions.

I. Assembly of stiffness matrix
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
- Objective: get the Mex function mx_assemble_sparse.mexa64
(Linux), mx_assemble_sparse.mexmaci64 (MacOS) or mx_assemble_sparse.mexw64
(Windows) (this function can be called directly in MATLAB)
- The codes are found in the directory "mx_FEM"
- The code providing the interface is mx_assemble_sparse.cpp 
- Require Armadillo C++ library (latest version tested with is 5.400.2, which
is included in the directory "Armadillo_packages")
- To compile the Mex function, one can either use the Makefile or the Matlab
script compile_assemble_mx.m
- Before using Makefile, ensure that: 
i) the Armadillo include file path is specified via the
compiler flag -I (see the CXXFLAGS variable in Makefile),
ii) the Armadillo library flag
-larmadillo and the gfortran library flag -lgfortran are specified,
iii) the mex compiler knows where the Armadillo libraries called 
librarmadillo* are. This can be achieved
by specifying -L<path to library> after the -larmadillo command. I usually add
that library path in the "mexopts.sh" file in the Matlab hidden directory
~/.matlab/<Matlab version>. In addition, one also need to provide that path to
the environment variable LD_LIBRARY_PATH (otherwise, one would receive a
runtime library error during the linking step),
iv) the "mexopts.sh" file is properly configured (this is mainly system
dependent but is generally more straightforward in Linux than MacOS)
- Type Make to compile the Mex function together with other C codes
- If you prefer to compile directly from inside MATLAB, you can use the
script compile_assemble_mx.m instead of the Makefile. Ensure that -I<path to
Armadillo include files> and -L<path to Armadillo libraries> are properly set.
As usual, -larmadillo, -lgfortran must also be specified. This method is
generally slower as all files are recompiled whether they have been modified or
not

II. Intersection finder between the channel described by NURBS and the mesh
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
- Objective: get the Mex functions 
			 (i) curve_edge_nurbs_nurbs_intersect.mex*
			 (ii) nurbsSurf_phy2param.mex*
- The codes are found in a separate directory called SISL
- Require installation of SISL library (source codes are in SISL-master
directory of the SISL directory)
- Can either use Makefile or compile_sisl_function.m
- Ensure that Mex knows where the header file sisl.h is
- Specify -lsisl when linking and ensure that the path to libsisl.a is known by
Mex

Common issues 
-------------
Matlab cannot find Armadillo shared library
- 	Either provide the path of the library to the environment variable LD_LIBRARY_PATH outside of MATLAB
	For example in the bashrc, add this line
		LD_LIBRARY_PATH = <path to armadillo shared library>:$LD_LIBRARY_PATH
- 	or to hard-code the runtime library path in the compiled binary, specify
	-Wl,-rpath,<path> in linking step

Directories/subdirectories
--------------------------
IGFEM-Curves-2D
- Armadillo_packages
- ChannelFiles
- M_channels
- M_error_analysis
- M_FEM
- M_geom_toolbox
- M_postprocessing
- M_preFEM
- M_tests
- mx_FEM
NURBS
- nurbs-toolbox
SISL
- SISL-masters

Support
-------
If you are having issues or would like to develop these codes further, please
contact the author via email.

Walkthrough 
----------- 
A few walkthrough examples are provided in the file User_manual.pdf in the
directory Manuals

Copyright
---------
2013-2017 University of Illinois, Urbana-Champaign. All rights reserved.

Author
------
Marcus Hwai Yik Tan
PhD, Theoretical and Applied Mechanics
