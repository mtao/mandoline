# Mandoline

A robust cut-cell mesh generator for arbitrary triangle-mesh inputs

## Dependencies
### Compilation
Mandoline depends on a C++17 enabled compiler
- [gcc](https://gcc.gnu.org) with version >= 8.2.0 should work.
- [cmake](https://cmake.org) with version >= 3.10.1 .
```bash
apt install git cmake build-essentials gcc
```

as well as a few other annoying things:

### Built-in Libraries
Mandoline depends on Michael Tao's messy ``[core](https://github.com/mtao/core)'' library, [libigl](https://github.com/libigl/libigl).

- core and libigl are included as submodules, and can be added by:
```bash
git submodule update --init --recursive
```

### Required External Libraries
- A relatively new verison of [Eigen3](https://eigen.tuxfamily.org) is required (3.3.7 or the development branch?). The version provided by libigl will not suffice.
- [Boost](https://boost.org) is required for threads.
- [Google Protobuf](https://developers.google.com/protocol-buffers/) is required for the default cutmesh serialization format.
- [Eigen](https://eigen.tuxfamily.org) is required for almost every bit of linear algebra. With visualization it currently needs the dev branch (requires [mercurial](https://www.mercurial-scm.org) and check ```build-utils/fetch_dependncies_ubuntu.sh``` for more details on how to build)

```bash
apt install libboost-thread-dev  libeigen3-dev protobuf-compiler 
```



In order to handle self-intersections Mandoline requires, on top of libigl,
- [CGAL](https://www.cgal.org)
- [gmp](https://gmplib.org)
- [mpfr](https://www.mpfr.org).

```bash
apt install git cmake libboost-thread-dev libmpfr-dev libmpfrc++-dev libcgal-dev libeigen3-dev protobuf-compiler mercurial debhelper libgl-dev libopenal-dev libglfw3-dev libsdl2-dev libbullet-dev libglm-dev
libmpfr-dev libmpfrc++-dev libcgal-dev
```

### Optional Libraries
- [OpenMP](https://www.openmp.org) is a really good idea.
- [Magnum](https://github.com/mosra) and [Corrade](https://github.com/mosra/corrade) for visualization.
If Magnum and Corrade are available, targets that depend on visualization can be enabled by ```-DUSE_OPENGL ON``` when running cmake.
Installing Magnum and Corrade can be a pain, as I use the imgui plugin and default installations dont have that available. If you want to use it you might want to use the following script with MANDOLINE_DIR with the directory mandoline is stored in:

```bash
apt install debhelper libgl-dev libopenal-dev libglfw3-dev libsdl2-dev libbullet-dev libglm-dev

mkdir mosra; pushd mosra;
for repo in corrade magnum magnum-integration;
do git clone https://github.com/mosra/$repo; done

pushd magnum-integration;
pwd
sed "s|IMDIR|${IMGUI_DIR}|g" ${MANDOLINE_DIR}/build-utils/magnum_integration.patch > magnum_integration_imgui.patch
patch -b ./package/debian/rules  magnum_integration_imgui.patch
popd #magnum integration

for repo in corrade magnum magnum-integration;
do pushd $repo;
    ln -s package/debian .
    # ignoring the dsc errors taht happen from magnum
    dpkg-buildpackage || true
    popd
done
popd
```
