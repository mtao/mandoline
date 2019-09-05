# Mandoline

A robust cut-cell mesh generator for arbitrary triangle-mesh inputs

## Dependencies
### Compilation
Mandoline depends on a C++17 enabled 

### Built-in Libraries
Mandoline depends on Michael Tao's messy ``[core](https://github.com/mtao/core)'' library, [libigl](https://github.com/libigl/libigl).

-core and libigl are included as submodules, and can be added by:
```bash
git submodule update --init --recursive
```

### Required External Libraries
- A relatively new verison of [Eigen3](https://eigen.tuxfamily.org) is required (3.3.7 or the development branch?). The version provided by libigl will not suffice.
- [Boost](https://boost.org) is required for threads.
- [OpenMP](https://www.openmp.org) is required.

In order to handle self-intersections Mandoline requires 
- [CGAL](https://www.cgal.org)
- [gmp](https://gmplib.org)
- [mpfr](https://www.mpfr.org).

### Optional Libraries
Magnum depends on [Magnum](https://github.com/mosra) and [Corrade](https://github.com/mosra/corrade) for visualization.
If Magnum and Corrade are available, targets that depend on visualization can be enabled by ```-DUSE_OPENGL ON``` when running cmake.
