# Mandoline
[![CircleCI](https://circleci.com/gh/mtao/mandoline/tree/master.svg?style=svg)](https://circleci.com/gh/mtao/mandoline/tree/master)
A robust cut-cell mesh generator for arbitrary triangle-mesh inputs


## Serialization
By default Mandoline uses [protobuf](https://developers.google.com/protocol-buffers/) for serialization, following the format of defined at [proto/cutmesh.proto](https://github.com/mtao/mandoline/blob/master/proto/cutmesh.proto). This is language agnostic and so in theory, so long as your language has protobuf support/bindings you can read the output of Mandoline. In fact a lot of initial validity testing on Mandoline were done in python.

## Examples/Tools
There are a few examples / tools that are provided.
They can be built by 
```bash
make example_name
```

### curve_drawing_2d
In ```examples/curve_drawing_2d``` there is an example of how to run the 2D version of Mandoline on a drawn curve.

After drawing a curve press `Make CCM` to build the cutmesh. By default it should show Random colors for each cut-cell, but that visualization mode can be changed by the Combo box at the top of the gui widget (regions = 1 color per region; Harmonic = the solution to a Poisson equation based off the boundary; Harmonic_RHS shows the right hands used in the Poisson equations).

### Boundary curve input/output
in `boundary_curve_to_cutmesh2` one can pass in a cutmesh file of the format defined in [2D Format](https://github.com/mtao/mandoline/wiki/2D-Format).


### make_cutmesh_gui
We provide a visual interface for creating cutmeshes from OBJ files whose code is found at  ```tools/make_cutmesh_gui.cpp```, which provides a gui widget for changing the grid resolution, grid bounding box, etc before making a cutmesh and saving it to disk. The obj file is a commandoline argument.
```bash
make_cutmesh_gui obj_filename.obj
```
### cutmesh_info
We provide a simple tool for inspecting high level details about a cutmesh file, whose code is found in ```tools/cutmesh_info.cpp```. The meat of it is found at ```src/tools/cutmesh_info.cpp``` for those curious.
```bash
cutmesh_info cutmesh_file.cutmesh
```


### Compilation
Mandoline depends on a C++17 enabled compiler
- [gcc](https://gcc.gnu.org) with version >= 8.2.0 or [clang](https://clang.llvm.org) with version >= 8 should work.
- [cmake](https://cmake.org) with version >= 3.11.1 .
```bash
apt install git cmake build-essential gcc
```

as well as a few other annoying things in the subsequent sections.

Once everything in subsequent sections has been installed, mandoline can be bulit by

```bash
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j5
```

### Built-in Libraries
Mandoline depends on a few libraries that it will fetch on its own:
- Michael Tao's messy ``[core](https://github.com/mtao/core)'' library,
- [libigl](https://github.com/libigl/libigl) for its ability to make intersection-free triangle meshes,


### Required External Libraries
- A relatively new verison of [Eigen3](https://eigen.tuxfamily.org) is required (3.3.7 or the development branch?). The version provided by libigl will not suffice. MAndoline might fetch its own if your Eigen is not new enough.
- [tbb](https://software.intel.com/en-us/tbb) for multithreading.
- [Google Protobuf](https://developers.google.com/protocol-buffers/) is required for the default cutmesh serialization format.
- [Eigen](https://eigen.tuxfamily.org) is required for almost every bit of linear algebra. If the currently available version is insufficeintly high mandoline's build will fetch a sufficiently new version

```bash
apt install libboost-thread-dev  libeigen3-dev protobuf-compiler 
```



In order to handle self-intersections Mandoline requires, on top of libigl,
- [CGAL](https://www.cgal.org)
- [gmp](https://gmplib.org)
- [mpfr](https://www.mpfr.org).

```bash
apt install libboost-thread-dev libmpfr-dev libmpfrc++-dev libcgal-dev 
libmpfr-dev libmpfrc++-dev libcgal-dev
```

### Optional Libraries
- [OpenMP](https://www.openmp.org) is a really good idea, Eigen uses it a [bit](https://eigen.tuxfamily.org/dox/TopicMultiThreading.html) but I'm trying to switch to `tbb` for most uses in Mandoline itself.
- [Magnum](https://github.com/mosra/magnum) and [Corrade](https://github.com/mosra/corrade) (version >= 2019.10) for visualization.
- [imgui](https://github.com/ocornut/imgui) for visualizing things as well.
Targets that depend on visualization can be enabled by ```-DUSE_OPENGL ON``` when running cmake. They will be automatically fetched.

OpenGL stuff probably requiers the following on ubuntu:
```bash
apt install libgl-dev libopenal-dev libglfw3-dev libsdl2-dev libbullet-dev libglm-dev
```

