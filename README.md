# Mandoline

A robust cut-cell mesh generator for arbitrary triangle-mesh inputs

## Dependencies
### Compilation
Mandoline depends on a C++17 enabled  compiler
- [gcc](https://gcc.gnu.org) with version >= 8.2.0 should work.
- [cmake](https://cmake.org) with version >= 3.10.1 .

### Built-in Libraries
Mandoline depends on Michael Tao's messy ``[core](https://github.com/mtao/core)'' library, [libigl](https://github.com/libigl/libigl).

-core and libigl are included as submodules, and can be added by:
```bash
git submodule update --init --recursive
```

### Required External Libraries
- A relatively new verison of [Eigen3](https://eigen.tuxfamily.org) is required (3.3.7 or the development branch?). The version provided by libigl will not suffice.
- [Boost](https://boost.org) is required for threads.
- [Google Protobuf](https://developers.google.com/protocol-buffers/) is required for the default cutmesh serialization format.

In order to handle self-intersections Mandoline requires, on top of libigl,
- [CGAL](https://www.cgal.org)
- [gmp](https://gmplib.org)
- [mpfr](https://www.mpfr.org).

### Optional Libraries
- [OpenMP](https://www.openmp.org) is a really good idea.
- [Magnum](https://github.com/mosra) and [Corrade](https://github.com/mosra/corrade) for visualization.
If Magnum and Corrade are available, targets that depend on visualization can be enabled by ```-DUSE_OPENGL ON``` when running cmake.
