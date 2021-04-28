# GeoLib - basic 2D and 3D geometry
[![pipeline status](https://git.irt-e.de/driverless/geolib/badges/master/pipeline.svg)](https://git.irt-e.de/driverless/geolib/-/commits/master)

This library contains basic vector arithmetic and helper functions. It is built to be extended with any other geometry based function used in the car/simulation.


## Dependencies

- Unix based OS


## How to install

Do not use the script given in the repo. It is used for CI/CD only.

The library will be installed in your global path, therefore you can
pick a directory outside your repository for the installation.
If you do so, add the library with a link to your dependency list in the
README.md. 

In your directory of choice:
```sh
git clone git@gitlab.com:irt-driverless/geolib.git
mkdir -p geolib/build && cd geolib/build
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=g++ ..
make
make install
```
If you want to uninstall, use `make uninstall` in the `build` directory.


## How to use in your project

Make sure you followed the install steps line by line.
In your CMakeLists.txt file, add the following lines after your target.
```cmake
# header files for all librarys
target_include_directories(<exe_target> PUBLIC /usr/local/include)

# the librarys you want to link
target_link_libraries(<exe_target> PUBLIC /usr/local/lib/libgeo.a)
```

## How to extend

If you have a geometry function that is needed for your package or node
1. clone this repo `git clone https://git.irt-e.de/driverless/geolib`
2. make a feature-branch `git checkout -b feature/<your_feature>`
3. add tests for your function (tests/tests.cpp)
4. implement your function
5. push your branch and make a merge-request for master

Your request will be checked and merged.

