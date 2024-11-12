# Linalg Library

Original library: [https://gitlab.com/enrico_facca/p1galerkin/](https://gitlab.com/enrico_facca/p1galerkin/)

* Clone the repository
```
git clone git@github.com:ncrescenzio/p1galerkin.git
```

* Build the library
```bash
cd globals && mkdir build
cmake -S . -B build
cmake --build build
```

* Install the library
```
cmake --install build --prefix /path/to/install/dir
```

## Dependencies

* [`globals` library](https://github.com/ncrescenzio/globals)
* [`linalg` library](https://github.com/ncrescenzio/linear_algebra)
* [`geometry` library](https://github.com/ncrescenzio/geometry)
* `lapack`
* `blas`
