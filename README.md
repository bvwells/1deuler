# 1deuler
[![Build Status](https://travis-ci.org/bvwells/1deuler.svg?branch=master)](https://travis-ci.org/bvwells/1deuler)

One dimensional Compressible Euler Equations solved with moving mesh approach described in the PhD thesis

*A moving mesh finite element method for the numerical solution of partial differential equations and systems.*

which can be found [here][1].

## Numerical Solution

The one-dimensional Euler Equations are solved using a moving mesh 
method which uses the monitor function ```M=u(x,t)``` in the moving mesh 
equations for the mesh velocity. The mesh is advanced forwards in time 
using a forward Euler time-stepping scheme. The solution to the Euler equations
are obtained by solving the Arbitrary Lagrangian Eulerian (ALE) formation
of the equations with a finite volume method. All the moving mesh equations are solved using
linear finite elements.

## Building and Developing

Developing locally requires Docker for Windows. Run the command

```
docker build -t 1deuler .
```

to build the docker image which pulls the gcc base image containing gfortran and maps the source code into the container.

Then run image with the command:

```
docker run -i -t -v /f/git/src/github.com/bvwells/1deuler:/app 1deuler
```

This command maps the local workspace into the running image so any changes made in the running image will be reflected on the local workspace.

Within the running image generate the make files for the release version by running the command:

```
cmake .
```

To build the debug version of the code run the command:

```
cmake -DCMAKE_BUILD_TYPE=Debug
```

Build the executable by running the command:

```
make
```

[1]: http://www.reading.ac.uk/nmsruntime/saveasdialog.aspx?lID=24080&sID=90294
