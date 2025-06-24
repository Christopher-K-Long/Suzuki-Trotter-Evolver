# Switching Function

So far we have considered the situation in which we know the control amplitudes but not the evolution of the state. The opposite problem is often of interest: which control amplitudes produce a desired evolution. This is the subject of optimal control. A common solution is to optimise for the control amplitudes with respect to a cost function. A building block of this optimisation procedure is the switching function: The variational derivative of the cost function with respect to the control amplitudes.

In what follows we will consider Mayer problems: problems in which the cost function $J$ is a function of the final state $\psi[\vec a(t);T]$. More specifically, we will take the cost function to be an expectation value of an observable $\hat O$ with respect to the final state:

$$
J[\vec a(t)]=J(\psi[\vec a(t);T])=\psi^\dagger[\vec a(t);T]\hat O\psi[\vec a(t);T].
$$(switching_function_1)

Using variational calculus we find the switching function is given by

$$
\phi_j(t)\coloneqq\frac{\delta J}{\delta a_j(t)} =2\operatorname{Im}\left(\psi^\dagger[\vec a(t);T] \hat OU(t\to T)H_j\psi[\vec a(t);t]\right),
$$(switching_function_2)

where $U(t\to T)$ is the unitary that evolves the system from time $t$ to time $T$.

However, our evolution is being approximated with a Suzuki-Trotter expansion—see [previous section](state_vector_evolution.md). Thus, in order to optimise our control amplitudes we should calculate our gradients by differentiating the Suzuki-Trotter expansion:

$$
\phi_j(n\Delta t)\approx\tilde \phi_j(n\Delta t)\coloneqq\frac{1}{\Delta t}\frac{\partial J}{\partial a_{nj}}
$$(switching_function_3)

$$
=\!2\operatorname{Im}\!\left(\psi^\dagger(T)
\hat O\!\!\left[\prod_{i>n}^N\prod_{k=1}^{\texttt{length}}
e^{-ia_{ik}H_k\Delta t}\right]\!\!\!
\left[\prod_{k=j}^{\texttt{length}}
e^{-ia_{nk}H_k\Delta t}\right]\!H_j\!\!
\left[\prod_{k=0}^{j-1}
e^{-ia_{nk}H_k\Delta t}\right]
\!\psi(\left[n-1\right]\Delta t)\right).
$$(switching_function_4)
where for numerical efficiency we replace $e^{-ia_{ik}H_k\Delta t}$ with $U_ke^{-ia_{ik}D_k\Delta t}U_k^\dagger$ as in [previous section](state_vector_evolution.md).

## Example

Numerically, we can compute $\tilde \phi_j(n\Delta t)$ using [``switching_function()``](../reference/structSuzuki__Trotter__Evolver_1_1UnitaryEvolver.rst#_CPPv4N22Suzuki_Trotter_Evolver14UnitaryEvolver18switching_functionE7DMatrixI7Dynamic6n_ctrlE7DMatrixI3dimXL1EEEd7DMatrixI3dim3dimE). As a worked example we will extend our Rabi oscillation example with a Pauli-y drive and we let

$$
\hat O=\begin{bmatrix}1&\phantom{-}1\\1&-1\end{bmatrix}.
$$

To simulate this we can execute the following program:

```{literalinclude} ../../examples/compute_switching_function.cpp
```

Output:

```
Expectation value: 0.577827
Switching function:
(-0.861075,0)   (2.42157,0)
(-0.981027,0)   (2.33944,0)
 (-1.09672,0)   (2.22863,0)
 (-1.20674,0)   (2.09098,0)
 (-1.30974,0)   (1.92885,0)
  (-1.4045,0)   (1.74503,0)
 (-1.48996,0)   (1.54262,0)
  (-1.5652,0)   (1.32497,0)
 (-1.62946,0)   (1.09558,0)
 (-1.68218,0)  (0.857956,0)
 (-1.72296,0)  (0.615554,0)
 (-1.75157,0)  (0.371674,0)
 (-1.76796,0)  (0.129392,0)
 (-1.77222,0) (-0.108504,0)
 (-1.76458,0)  (-0.33956,0)
  (-1.7454,0) (-0.561681,0)
 (-1.71515,0)  (-0.77315,0)
 (-1.67436,0) (-0.972622,0)
 (-1.62366,0)  (-1.15911,0)
  (-1.5637,0)  (-1.33197,0)
```

The switching function output is a tuple. The first entry is the cost function (expectation value). Note we can compute the expectation value alone using [``evolved_expectation_value()``](../reference/structSuzuki__Trotter__Evolver_1_1UnitaryEvolver.rst#_CPPv4N22Suzuki_Trotter_Evolver14UnitaryEvolver25evolved_expectation_valueE7DMatrixI7Dynamic6n_ctrlE7DMatrixI3dimXL1EEEd7DMatrixI3dim3dimE). The second entry in the tuple is the switching function itself as a $\texttt{length}\times N$ matrix.

---
[Previous](state_vector_evolution.md) | [Next](examples.md)