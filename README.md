# Hobak - A Library for Squashing Things
![HOBAK Logo](./data/HOBAK_logo.png)

This is the ***Hobak*** library, which accompanies the SIGGRAPH 2022 Course *[Dynamic Deformables: Implementation and Production Practicalities (Now With Code!)](http://www.tkim.graphics/DYNAMIC_DEFORMABLES/)* 
by Theodore Kim and David Eberle. This library is intended to illustrate many of the concepts from the course as clearly as possible,
and as a consequence is currently **not optimized for performance**. Readers beware: anybody claiming a 5X speed win over this
code ain't actually claiming that much.

***Hobak*** is Korean for *squash* or *pumpkin*. This is library for squashing things.

## Quick Start

***Hobak*** is designed to build and run with the smallest number of dependencies possible. All the necessary external libraries
are already included in the `HOBAK/ext` directory. Hopefully, you will not need to download anything else.

The fastest way to get up and running is to pick a build architecture, and call `make` on it. From the top level `HOBAK` directory,
to build on Linux, just call `make linux`. It should put everything where it needs to go, build the `simulateScenes` project, 
and run a *Bunny Drop* simulation. 

If you want to build on a Mac, just call `make mac`. This one might give you a build error saying it can't find GLUT, in which 
case you need to install [XQuartz](https://www.xquartz.org/), because apparently Mac is just toooo goooood to have GLUT installed 
by default.

If you have a version of `g++` installed that supports OpenMP (the default Mac compiler, `clang`, does not), you can alternatively 
build everything on the Mac by calling `make mac_omp`. It will expect that `g++` already points to an OpenMP-capable compiler.

## Building Other Projects

To build a project in ***Hobak***, you first need to pick an architecture. If you called `make linux` or `make mac`, then it
already picked for you.

If you want to pick by hand, go to `HOBAK/projects/`, and copy one of the `include_top` files to `include_top.mk`. 
For example, to build on Linux, call `cp include_top.linux include_top.mk`. Again, if you already ran `make linux`, this was 
already done for you.

From there, if you want to build a specific project, go to its directory and call `make depend`, followed by `make`.
For example, to build the regression test suite, do:

    cd projects/regressionTests
    make depend
    make

You only need to call `make depend` once. It will build the dependency files, i.e. which files depend on which headers.
After calling it once, you should only need to call `make` from then on. You should only need to run `make depend` again if 
you change the dependencies, such as adding a new file to the `Makefile`, or changing a header.

## Running

Everything should be run from the `./bin` directory. So, to build and run `regressionTests`, do:

    cd projects/regressionTests
    make depend
    make
    cd ../../bin
    ./regressionTests

Projects that pop up an OpenGL window, like `simulateScenes` or `replayScenes` support the following keyboard commands:

 - `a` will start and stop the animation.
 - `q` will quit, and write out a JSON file of simulation data and a QuickTime movie of the animation. You need to have [FFMPEG](https://ffmpeg.org/)
installed for the movie write to succeed.
 - `Q` or `Esc` will quit without writing anything out.

## Unit and Regression Testing

Some basic unit and regression testing has been implemented using the *[Catch2](https://github.com/catchorg/Catch2)* library. The following
should build and run the units:

    cd projects/unitTests
    make depend
    make
    cd ../../bin
    ./unitTests

## What's Missing

A bunch of features are missing, mostly because my coding time is limited. If you want to see these features
implemented, either add them yourself or give me a grant.

 - **Continuous Collision Detection (CCD)** and **Global Intersection Analysis (GIA)**: If something weird happens
 under violent collisions, e.g. triangles get snagged, this is probably why.

 - **Line Search**: If the simulation explodes, this is probably why.

 - **Shells and Strands**: For now, you can fake it with really thin tet meshes.

## External Libraries

***Hobak*** uses the following external libraries, all of which are included in the download, so that we don't run into versioning
issues a year from now.

 - [Eigen](https://eigen.tuxfamily.org)
 - [Catch2](https://github.com/catchorg/Catch2)
 - [RapidJSON](https://rapidjson.org)
 - [Spectra](https://spectralib.org)
 - [GLVU](http://www.cs.unc.edu/~walk/software/glvu/)

<!---
## The Tobj file format

Several tet mesh files are provided for you in the *./data/* directory, such as a simple cube at different resolutions, 
and the Bunny at several resolutions. These are stored in a flat, text-based *.tobj* format, which is the same as the 
popular *[.obj](https://en.wikipedia.org/wiki/Wavefront_.obj_file)* format.

Vertices are specified the same way as before, and are numbered according to the order that they appear in the file

    v 1.0 2.0 1.0

The only new addition is a tetrahedron tag has also been added, denoted *t*, which indicates with vertices make up
a tetrahedron

    t 100 341 82 900

I sure hope the *t* tag wasn't already taken in the .obj specification. There wasn't any conflict that I could see.

To make your life slightly easier, I have also included a [Gmsh](http://gmsh.info) to Tobj converter, under the project
*gmshToTobj*. It uses the Gmsh loader from Qingnan Zhou's [PyMesh](https://github.com/PyMesh/PyMesh) library.
--->
