# Suzuki-Trotter-Evolver

A C++ header-only library for evolving states under the Schrödinger equation using first-order Suzuki-Trotter and computing switching functions.

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17116329.svg)](https://doi.org/10.5281/zenodo.17116329) [![SWH](https://archive.softwareheritage.org/badge/origin/https://doi.org/10.5281/zenodo.17116329/)](https://archive.softwareheritage.org/browse/origin/?origin_url=https://doi.org/10.5281/zenodo.17116329)

[![Build and Test Ubuntu-latest](https://github.com/Christopher-K-Long/Suzuki-Trotter-Evolver/actions/workflows/ubuntu_latest.yml/badge.svg)](https://github.com/Christopher-K-Long/Suzuki-Trotter-Evolver/actions/workflows/ubuntu_latest.yml) [![Build and Test macOS-latest](https://github.com/Christopher-K-Long/Suzuki-Trotter-Evolver/actions/workflows/macos-latest.yml/badge.svg)](https://github.com/Christopher-K-Long/Suzuki-Trotter-Evolver/actions/workflows/macos-latest.yml) [![Build and Test Windows-latest](https://github.com/Christopher-K-Long/Suzuki-Trotter-Evolver/actions/workflows/windows-latest.yml/badge.svg)](https://github.com/Christopher-K-Long/Suzuki-Trotter-Evolver/actions/workflows/windows-latest.yml)

## Installation

Suzuki-Trotter-Evolver can be installed as follows:

```bash
git clone https://github.com/Christopher-K-Long/Suzuki-Trotter-Evolver
cd Suzuki-Trotter-Evolver
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release # set optional flags here
sudo cmake --build build --config Release --target install -j $(nproc)
```

There are several optional flags that can be set:
- ``-DSuzuki-Trotter-Evolver_INSTALL_LIBRARY=ON/OFF``(Suzuki-Trotter-Evolver is installed if set to ``ON``)
- ``-DSuzuki-Trotter-Evolver_BUILD_TESTING=ON/OFF`` (The unit tests and examples are build and run if set to ``ON``)
- ``-DSuzuki-Trotter-Evolver_BUILD_DOCS=ON/OFF`` (The documentation is build if set to ``ON``)

These optional flags should be appended to ``cmake -S ../.. -B build``. All three default to ``ON`` if the project is top level (i.e. if the project is not being build as a dependency of another project).

### Requirements

Runs on Linux, macOS, and Windows. Requires:

- [Eigen3](https://eigen.tuxfamily.org/)

#### Additional requirements for testing

- [Catch2](https://github.com/catchorg/Catch2)

#### Additional requirements for building documentation

- [doxygen](https://doxygen.nl/)
- [Sphinx](https://www.sphinx-doc.org/)
- [Furo](https://github.com/pradyunsg/furo)
- [Exhale](https://exhale.readthedocs.io)
- [Breathe](https://breathe.readthedocs.io)
- [MyST-Parser](https://myst-parser.readthedocs.io/)

## Documentation

Documentation, including worked examples, can be found at: [https://Suzuki-Trotter-Evolver.readthedocs.io](https://Suzuki-Trotter-Evolver.readthedocs.io)

## Source Code

Source code can be found at: [https://github.com/Christopher-K-Long/Suzuki-Trotter-Evolver](https://github.com/Christopher-K-Long/Suzuki-Trotter-Evolver)

A mirror can be found at: [https://gitlab.com/Christopher-K-Long/Suzuki-Trotter-Evolver](https://gitlab.com/Christopher-K-Long/Suzuki-Trotter-Evolver)

Please submit all [pull requests](https://github.com/Christopher-K-Long/Suzuki-Trotter-Evolver/pulls), [issues](https://github.com/Christopher-K-Long/Suzuki-Trotter-Evolver/issues), [discussions](https://github.com/Christopher-K-Long/Suzuki-Trotter-Evolver/discussions), and [vulnerability reports](https://github.com/Christopher-K-Long/Suzuki-Trotter-Evolver/security) to the [GitHub](https://github.com/Christopher-K-Long/Suzuki-Trotter-Evolver) repository.

Releases from this repository are assigned DOIs and  can be found at [https://doi.org/10.5281/zenodo.17116329](https://doi.org/10.5281/zenodo.17116329). The DOI for all releases is [10.5281/zenodo.17116329](https://doi.org/10.5281/zenodo.17116329). Additionally, the releases are archived to [https://archive.softwareheritage.org/swh:1:dir:364f1095abecf6ebf480fbccacfe23f0a10dce63](https://archive.softwareheritage.org/swh:1:dir:364f1095abecf6ebf480fbccacfe23f0a10dce63).

## Version and Changes

The current version is [`1.1.0`](ChangeLog.md#release-110). Please see the [Change Log](ChangeLog.md) for more details. Suzuki-Trotter-Evolver uses [semantic versioning](https://semver.org/).
