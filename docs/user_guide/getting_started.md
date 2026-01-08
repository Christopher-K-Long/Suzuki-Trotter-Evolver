# Getting Started

## What is Suzuki-Trotter-Evolver

Suzuki-Trotter-Evolver is a C++ header-only library for evolving states under the Schrödinger equation,

$$
\dot\psi=-iH\psi,
$$(getting_started_1)

using the first-order Suzuki-Trotter expansion and computing switching functions. This library assumes Hamiltonians take the form:

$$
H(t)=H_0+\sum_{j=1}^{\textrm{length}}a_j(t)H_j,
$$(getting_started_2)

where $H_0$ is called the drift Hamiltonian and $\left\{H_j\right\}_{j=1}^{\textrm{length}}$ are the control Hamiltonians which are modulated by the control amplitudes $\left\{a_j(t)\right\}_{j=1}^{\textrm{length}}$.
Suzuki-Trotter-Evolver integrates the Schrödinger equation by diagonalising the drift and control
Hamiltonians ahead of time so that no matrix exponentials are required in the Suzuki-Trotter expansion. This
yields an upfront cost of $O(\dim^3)$ where $\dim$ is the dimension of the vector space the Hamiltonians act upon. However, integrating the Schrödinger equation will scale with $\dim$ as $O(\dim^2)$ opposed to $O(\dim^3)$.

### Is Suzuki-Trotter-Evolver for you?
As we have just discussed has a large upfront cost in order to reduce the integration cost. Thus, if you want to perform many integrations of the Schrödinger equation for fixed drift and control Hamiltonians (with potential different control amplitudes) then Suzuki-Trotter-Evolver is for you. However, if you want to perform only a few integrations or even a single integration then the upfront cost may outweigh the better scaling for each integration.

## Installation

Suzuki-Trotter-Evolver can be installed as follows:

```bash
git clone https://github.com/Christopher-K-Long/Suzuki-Trotter-Evolver
cd Suzuki-Trotter-Evolver
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release # set optional flags here
sudo cmake --build build --config Release --target install -j $(nproc)
```

There are several optional flags that can be set:
- ``-Suzuki-Trotter-Evolver_INSTALL_LIBRARY=ON/OFF`` (Suzuki-Trotter-Evolver is installed if set to ``ON``)
- ``-DSuzuki-Trotter-Evolver_BUILD_TESTING=ON/OFF`` (The unit tests and examples are build and run if set to ``ON``)
- ``-DSuzuki-Trotter-Evolver_BUILD_DOCS=ON/OFF`` (The documentation is build if set to ``ON``)

These optional flags should be appended to ``cmake -S ../.. -B build``. All three default to ``ON`` if the project is top level (i.e. if the project is not being build as a dependency of another project).


## Quick Start

A good starting point in experimenting with a new integrator for quantum systems is to try and reproduce a Rabi Oscillation. Consider the Hamiltonian

$$
H(t)=\frac{1}{2}vZ+f\cos(vt)X,
$$(getting_started_3)

where $X$ and $Z$ are the Pauli-$x$ and -$z$ matrices, respectively. This Hamiltonian produces on resonance Rabi-oscillations at a frequency $f$. Suppose we initialise a in the state $\psi(0)=\begin{pmatrix}1&0\end{pmatrix}^\intercal$ and want to know the state $\psi(t=T)$ at a time $t=T$. The following program will print this state to the console:

```{literalinclude} ../../examples/Rabi_oscillation.cpp
```

---
[Previous](overview.md) | [Next](state_vector_evolution.md)