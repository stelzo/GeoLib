# GeoLib - basic 2D and 3D geometry

This library contains basic vector arithmetic and helper functions. It is build to be extended with any other geometry based function used in the car/simulation.

## How to use

Do not use the script given in the repo. It is used for CI/CD only.

```sh
git clone git@git.irt-e.de:driverless/geolib.git
mkdir -p geolib/build && cd geolib/build
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=g++ ..
make
make install
```

If you want to uninstall, use `make uninstall` in the `build` directory.


## How to extend

If you have a geometry function that is needed for your package or node
1. clone this repo `git clone git@git.irt-e.de:driverless/geolib.git`
2. make a feature-branch `git checkout -b feature/<your_feature>`
3. add tests for your function (tests/tests.cpp)
4. implement your function
5. push your branch and make a merge-request for master

Your request will be checked and merged.

