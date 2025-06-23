# State Vector Evolution

We saw in [Quick Start](getting_started.md#quick-start) how we could implement a simulation of a Rabi Oscillation. Here we will first break down the error in the simulation and how to control this error before moving onto other ways we can propagate state vectors.

## Integration method and errors

As the name suggests this library uses the Suzuki-Trotter expansion to integrate the Schrödinger equation,

$$
\dot U=-iHU,
$$(state_vector_evolution_1)

with Hamiltonians of the form

$$
H(t)=H_0+\sum_{j=1}^{\textrm{length}}a_j(t)H_j.
$$(state_vector_evolution_2)

The first step we take to integrating the Schrödinger equation is approximating the Hamiltonian as constant over a period of time $\Delta t$. The error induced by this approximation can be derived as follows: First, integrate Eq. {eq}`state_vector_evolution_1`,

$$
U(t)=I-i\int_0^t\textrm{d}t_1H(t_1)U(t_1).
$$(state_vector_evolution_3)

Now recursively substituting the left-hand side of Eq. {eq}`state_vector_evolution_3` into the right-hand-side we find

$$
U(\Delta t)=I-i\int_0^{\Delta t}\textrm{d}t_1H(t_1)-\int_0^{\Delta t}\textrm{d}t_1\int_0^{t_1}\textrm{d}t_2H(t_1)H(t_2)+\mathcal O(\Delta t^3),
$$(state_vector_evolution_4)

where we have only retained terms up to order $\Delta t^2$. If we instead treated the Hamiltonian as constant over the time $\Delta t$ we would find

$$
U'(\Delta t)=e^{-iH(0)\Delta t}=I-iH(0)\Delta t-\frac{1}{2}H^2(0)\Delta t^2+\mathcal O(\Delta t^3).
$$(state_vector_evolution_5)

By Maclaurin expanding $H(t)$, we find the additive error induced by this approximation as

$$
U'(\Delta t)-U(\Delta t)=\frac{1}{2}i\dot H(0)\Delta t^2+\mathcal O(\Delta t^3).
$$(state_vector_evolution_6)

Thus, we can integrate the Schrödinger equation for a time $N\Delta t$ as:

$$
U(N\Delta T)=\prod_{n=0}^{N-1}e^{-iH(n\Delta t)\Delta t}+\mathcal O\left(\Delta t^2\sum_{n=0}^{N-1}\left\lvert\left\lvert\dot H(n\Delta t)\right\rvert\right\rvert\right).
$$(state_vector_evolution_7)

Next we use the first-order Suzuki-Trotter expansion to expand each exponential:[^1]

$$
e^{-iH(n\Delta t)\Delta t}=\prod_{j=0}^{\textrm{length}}e^{-ia_{nj}H_j\Delta t}+\mathcal O\left(\Delta t^2\sum_{j,k=0}^{\textrm{length}}a_{nj}a_{nk}\left\lvert\left\lvert\left[H_j,H_k\right]\right\rvert\right\rvert\right)
$$(state_vector_evolution_8)

where $a_{nj}=a_j(n\Delta t)$ and for notational ease we set $a_{n0}=1$.

Combining these results we find:

$$
U(N\Delta t)=\prod_{i=1}^N\prod_{j=0}^{\textrm{length}}e^{-ia_{ij}H_j\Delta t}+\mathcal E
$$(state_vector_evolution_9)

where

$$
\mathcal E=\mathcal O\left(
    \Delta t^2\left[\sum_{i=1}^N\sum_{j=1}^{\textrm{length}}\dot a_{ij}\left\lvert\left\lvert H_j\right\rvert\right\rvert+\sum_{i=1}^N\sum_{j,k=0}^{\textrm{length}}a_{ij}a_{ik}\left\lvert\left\lvert[H_j,H_k]\right\rvert\right\rvert\right]
    \right),
$$(state_vector_evolution_10)

and $\dot a_{nj}\coloneqq\dot a_j(n\Delta t)$. By making the following definitions

$$
\omega\coloneqq\max_{\substack{i\in\left[1,N\right]\\j\in\left[1,\textrm{length}\right]}}\left|\dot a_{ij}\right|,
$$(state_vector_evolution_11)

$$
\alpha\coloneqq\max_{\substack{i\in\left[1,N\right]\\j\in\left[0,\textrm{length}\right]}}\left|a_{ij}\right|,
$$(state_vector_evolution_12)

$$
E\coloneqq\max_{j\in\left[0,\textrm{length}\right]}\left\lvert\left\lvert H_j\right\rvert\right\rvert,
$$(state_vector_evolution_13)

we can simplify the error term to

$$
\mathcal E=\mathcal O\left(N\Delta t^2\textrm{length}\left[\omega E+\alpha^2+E^2\right]\right).
$$(state_vector_evolution_14)

Note the error is quadratic in $\Delta t$ but linear in $N$. We can also view this as being linear in $\Delta t$ and linear in total evolution time $N\Delta t$. Additionally, by Nyquist's theorem this asymptotic error scaling will not be achieved until the time step $\Delta t$ is smaller than $\frac{1}{2\Omega}$ where $\Omega$ is the largest energy or frequency in the system. In our [Rabi oscillation example](getting_started.md#quick-start) $\Omega=\max\left\{v,f\right\}$.

Finally, we can diagonalise each term in the Hamiltonian, $H_j=U_jD_jU_j^\dagger$ where $U_j$ is a unitary and $D_j$ is diagonal. Substituting this diagonalised form into Eq. {eq}`state_vector_evolution_9` we find

$$
U(N\Delta t)=\prod_{n=1}^N\prod_{j=0}^{\textrm{length}}U_je^{-ia_{nj}D_j\Delta t}U_j^\dagger+\mathcal E.
$$(state_vector_evolution_15)

When we act Eq. {eq}`state_vector_evolution_15` on a state vector $\psi(0)$,

$$
\psi(N\Delta t)=\prod_{n=1}^N\prod_{j=0}^{\textrm{length}}U_je^{-ia_{nj}D_j\Delta t}U_j^\dagger\psi(0)+\mathcal E,
$$(state_vector_evolution_16)

we can see how the acceleration is achieved by diagonalising the Hamiltonians. The matrix exponential $e^{-ia_{nj}H_j\Delta t}$ takes $O(\dim^3)$ time where $\dim$ is the dimension of the vector space acted upon by the Hamiltonians. Whereas, exponentiating the diagonal matrix $-ia_{nj}D_j\Delta t$ takes $O(\dim)$ time. The runtime is now dominated by the $O(\dim^2)$ matrix vector multiplications.

### Summary

| Property               | Scaling                                                                               |
| ---------------------- | ------------------------------------------------------------------------------------- |
| Initialisation runtime | $O(\textrm{length}\times\dim^3)$                                                      |
| Integrator runtime     | $O(N\times\textrm{length}\times\dim^2)$                                               |
| Integrator error       | $\mathcal O\left(N\Delta t^2\textrm{length}\left[\omega E+\alpha^2+E^2\right]\right)$ |

## Other propagation methods

In [Quick Start](getting_started.md#quick-start) we integrated a Rabi oscillation using [``propagate()``](../reference/structSuzuki__Trotter__Evolver_1_1UnitaryEvolver.rst#_CPPv4N22Suzuki_Trotter_Evolver14UnitaryEvolver9propagateE7DMatrixI7Dynamic6n_ctrlE7DMatrixI3dimXL1EEEd). Alternatively, we could use [``propagate_all()``](../reference/structSuzuki__Trotter__Evolver_1_1UnitaryEvolver.rst#_CPPv4N22Suzuki_Trotter_Evolver14UnitaryEvolver13propagate_allE7DMatrixI7Dynamic6n_ctrlE7DMatrixI3dimXL1EEEd) to collect the whole evolution of the state, $\left(\psi(n\Delta t)\right)_{n=0}^N$, simply by replacing [``propagate()``](../reference/structSuzuki__Trotter__Evolver_1_1UnitaryEvolver.rst#_CPPv4N22Suzuki_Trotter_Evolver14UnitaryEvolver9propagateE7DMatrixI7Dynamic6n_ctrlE7DMatrixI3dimXL1EEEd) with [``propagate_all()``](../reference/structSuzuki__Trotter__Evolver_1_1UnitaryEvolver.rst#_CPPv4N22Suzuki_Trotter_Evolver14UnitaryEvolver13propagate_allE7DMatrixI7Dynamic6n_ctrlE7DMatrixI3dimXL1EEEd).

Further, if we wish to propagate multiple states one solution is to wrap [``propagate()``](../reference/structSuzuki__Trotter__Evolver_1_1UnitaryEvolver.rst#_CPPv4N22Suzuki_Trotter_Evolver14UnitaryEvolver9propagateE7DMatrixI7Dynamic6n_ctrlE7DMatrixI3dimXL1EEEd) in a for loop:

```cpp
// Each column is an initial state
MatrixXcd initial_states {{1, 0, 1/std::sqrt(2)},
                          {0, 1, 1/std::sqrt(2)}};

// Prepare the matrix to store the output
MatrixXcd output_states(2, initial_states.cols());

// Propagate each input state
for (int j = 0; j < initial_states.cols(); j++) {
    output_states.col(j) = evolver.propagate(ctrl_amp,
                                             initial_states.col(j),
                                             dt);
}
```

However, a more concise and efficient code uses [``propagate_collection()``](../reference/structSuzuki__Trotter__Evolver_1_1UnitaryEvolver.rst#_CPPv4I_iEN22Suzuki_Trotter_Evolver14UnitaryEvolver20propagate_collectionE7DMatrixI3dim1lE7DMatrixI7Dynamic6n_ctrlE7DMatrixI3dim1lEd):

```cpp
// Each column is an initial state
MatrixXcd initial_states {{1, 0, 1/std::sqrt(2)},
                          {0, 1, 1/std::sqrt(2)}};

// Propagate each input state
MatrixXcd output_states = evolver.propagate_collection(ctrl_amp,
                                                       initial_states,
                                                       dt);
```

---
[Previous](getting_started.md) | [Next](switching_function.md)


## References

[^1]: Childs, A. M., Su, Y., Tran, M. C., Wiebe, N., & Zhu, S. (2021). Theory of Trotter Error with Commutator Scaling. Phys. Rev. X, 11, 011020. [doi:10.1103/PhysRevX.11.011020](https://www.doi.org/10.1103/PhysRevX.11.011020)