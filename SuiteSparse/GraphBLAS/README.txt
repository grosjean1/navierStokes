SuiteSparse:GraphBLAS, Timothy A. Davis, (c) 2017, All Rights Reserved.
http://suitesparse.com   See GraphBLAS/Doc/License.txt for license.

VERSION 1.1.0, Dec 1, 2017

SuiteSparse:GraphBLAS is an full implementation of the GraphBLAS standard,
which defines a set of sparse matrix operations on an extended algebra of
semirings using an almost unlimited variety of operators and types.  When
applied to sparse adjacency matrices, these algebraic operations are equivalent
to computations on graphs.  GraphBLAS provides a powerful and expressive
framework for creating graph algorithms based on the elegant mathematics of
sparse matrix operations on a semiring.

See https://graphblas.org for more information on GraphBLAS, including the
GraphBLAS C API.  See the user guide in the Doc/ folder for documentation on
the SuiteSparse implementation of GraphBLAS.

QUICK START: To compile, run several demos, and install, do these commands in
this directory:

    make
    sudo make install

Please be patient; some files can take several minutes to compile.  Requires an
ANSI C11 compiler, so cmake will fail if your compiler is not C11 compliant.
See the User Guide in Doc/*pdf for directions on how to use another compiler.

The output of the demo programs will be compared with their expected output.

To remove all compiled files:

    make clean

To compile the library without running the demos or installing it:

    make library

NOTE: this package has not yet been ported to Windows.  It uses cmake to build
the package so porting to Windows should be straight-forward (in progress).

--------------------------------------------------------------------------------
Files and folders in this GraphBLAS directory:

CMakeLists.txt  cmake instructions to compile GraphBLAS

Makefile        a very simple Makefile that relies on cmake

Demo            a set of demos on how to use GraphBLAS

Doc             SuiteSparse:GraphBLAS User Guide and license

Include         user-accessible include file, GraphBLAS.h

install         default installation folder for compiled GraphBLAS library
                (created by "make install", does not appear in the distribution)

Lib             for the SuiteSparse:GraphBLAS compiled library

Test            Extensive tests, not meant for general usage.  To compile
                SuiteSparse:GraphBLAS and test in MATLAB, go to this directory
                and type gbmake;testall in MATLAB.

Makefile        to compile the SuiteSparse:GraphBLAS library and demos

README.txt      this file

Source          source files of the SuiteSparse:GraphBLAS library.

Tcov            test coverage, requires MATLAB

build           build directory, intially empty

--------------------------------------------------------------------------------

SPEC: This version fully conforms to GraphBLAS C API Specification 1.1.0.
It includes several additional functions and features as extensions to the
spec.  These extensions are tagged with the keyword SPEC: in the code and in
the User Guide, and in the Include/GraphBLAS.h file.

