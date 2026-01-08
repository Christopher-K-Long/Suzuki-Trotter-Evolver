# Running Tests

You can run the unit tests yourself during the [build process](getting_started.md#installation) by setting the flag ``Suzuki-Trotter-Evolver_BUILD_TESTING`` to ``ON`` as follows:

```bash
git clone https://github.com/Christopher-K-Long/Suzuki-Trotter-Evolver
cd Suzuki-Trotter-Evolver
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DSuzuki-Trotter-Evolver_BUILD_TESTING=ON
cmake --build build -config Release --target install -j $(nproc)
```

The unit tests and the [examples](examples.md) will be run after they are built. If you wish to run the unit tests again this can be done as follows:

```bash
cd build
make test
```

or by executing

```bash
cd build
ctest
```

[Previous](examples.md)