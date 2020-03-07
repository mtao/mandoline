# Mandoline
[![CircleCI](https://circleci.com/gh/mtao/mandoline/tree/master.svg?style=svg)](https://circleci.com/gh/mtao/mandoline/tree/master)
A robust cut-cell mesh generator for arbitrary triangle-mesh inputs

## Dependencies


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
make make_cutmesh_gui -j5
```
where ```make_cutmesh_gui``` is a little UI for making cutmesehs from triangle meshes.

### Built-in Libraries
Mandoline depends on a few libraries that it will fetch on its own:
- Michael Tao's messy ``[core](https://github.com/mtao/core)'' library,
- [libigl](https://github.com/libigl/libigl) for its ability to make intersection-free triangle meshes,


### Required External Libraries
- A relatively new verison of [Eigen3](https://eigen.tuxfamily.org) is required (3.3.7 or the development branch?). The version provided by libigl will not suffice. MAndoline might fetch its own if your Eigen is not new enough.
- [Boost](https://boost.org) is required for threads.
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
- [OpenMP](https://www.openmp.org) is a really good idea.
- [Magnum](https://github.com/mosra/magnum) and [Corrade](https://github.com/mosra/corrade) (version >= 2019.10) for visualization.
- [imgui](https://github.com/ocornut/imgui) for visualizing things as well.
Targets that depend on visualization can be enabled by ```-DUSE_OPENGL ON``` when running cmake. They will be automatically fetched

OpenGL stuff probably requiers the following:
```bash
apt install libgl-dev libopenal-dev libglfw3-dev libsdl2-dev libbullet-dev libglm-dev
```

