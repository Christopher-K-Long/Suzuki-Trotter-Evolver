/**
    \file
    \brief An integrator for first-order homogeneous ordinary differential
    equations that assumes the solution is unitary and associated types and
    functions.
*/

#include <complex>
#include <cmath>
#include <vector>
#include <Eigen/Eigenvalues>
#include <Eigen/SparseCore>
using std::complex;
using std::vector;
using Eigen::Dynamic;
using std::pow;
using std::conj;

/**
    A C++ header-only library for evolving states under the Schrödinger equation
    using first-order Suzuki-Trotter and computing switching functions.
 */
namespace Suzuki_Trotter_Evolver {

/**
    A dense `n`x`m` matrix of complex doubles.
    @tparam n The number of rows.
    @tparam m The number of columns.

    Also see \ref SMatrix.
*/
template<int n = Dynamic, int m = Dynamic> using DMatrix =
    Eigen::Matrix<complex<double>, n, m>;

/**
    A sparse matrix of complex doubles.

    Also see \ref DMatrix.
*/
typedef Eigen::SparseMatrix<complex<double>> SMatrix;

/**
    Computes the gate infidelity between a gate and a target gate:
    @f[
    \mathcal I(\texttt{gate}, \texttt{target})
    \coloneqq 1-\frac{
    \left|\Tr\left[\texttt{target}^\dagger \cdot \texttt{gate}\right]\right|^2
    +\texttt{dim}}{\texttt{dim}(\texttt{dim}+1)},
    @f]
    where \f$\texttt{dim}\f$ is the dimension of the Hilbert space the gates act
    upon.

    @param gate The gate to compute the infidelity of.
    @param target The target gate to compute the infidelity with respect to.
    @return The gate infidelity:
    \f$\mathcal I(\texttt{gate}, \texttt{target})\f$.
*/
template <typename Derived>
double unitary_gate_infidelity(const Eigen::MatrixBase<Derived>& gate,
                               const Eigen::MatrixBase<Derived>& target) {
    Eigen::Index d = target.rows();
    return 1-(std::pow(std::abs((target.adjoint()*gate).trace()),2)+d) /
             (d*(d+1));
};

/**
    A struct to store the diagonalised drift and control Hamiltonians. On
    initialisation the Hamiltonians are diagonalised and the eigenvectors and
    values stored. This initial diagonalisation may be slow and takes
    \f$O(\texttt{dim}^3)\f$
    time for a
    \f$\texttt{dim}\times \texttt{dim}\f$
    Hamiltonian. However, it  allows each step of the Suzuki-Trotter expansion
    to be implimented in
    \f$O(\texttt{dim}^2)\f$
    time with matrix multiplication and only scalar exponentiation opposed to
    matrix exponentiation which takes
    \f$O(\texttt{dim}^3)\f$
    time.

    @tparam n_ctrl The number of control Hamiltonians.
    @tparam dim The dimension of the vector space the Hamiltonians act upon.
    @tparam Matrix The type of matrix to use. `Matrix` must take the value
    \ref DMatrix<dim, dim> or \ref SMatrix for dense or sparse matrices,
    respectively.
*/
template<int n_ctrl = Dynamic,
         int dim = Dynamic,
         typename Matrix = DMatrix<dim, dim>>
struct UnitaryEvolver {
    private: const complex<double> minus_i = complex<double>(0.0, -1.0);
    public:

    /**
        The dimension of rows in each control Hamiltonian multiplied by the
        number of control Hamiltonians. This is the number of rows in the
        ``control_hamiltonians`` argument of
        ``UnitaryEvolver::UnitaryEvolver()``.
    */
    static const int dim_x_n_ctrl = (dim == Dynamic || n_ctrl == Dynamic)
                                    ? Dynamic
                                    : dim * n_ctrl;
    /**
        The number of control Hamiltonians
    */
    size_t length;

    /**
        The eigenvalues,
        \f$\operatorname{diag}(D_0)\f$,
        of the drift Hamiltonian:
        \f$H_0=U_0D_0U_0^\dagger\f$.
    */
    Eigen::Array<complex<double>, dim, 1> d0;

    /**
        The eigenvalues,
        \f$\left(\operatorname{diag}(D_i)\right)_{i=1}^{\texttt{length}}\f$,
        of the control Hamiltonians:
        \f$H_i=U_iD_iU_i^\dagger\f$
        for all
        \f$i\in\left[\texttt{length}\right]\f$.

        Also see \ref ds, \ref u0, and \ref u0_inverse.
    */
    vector<Eigen::Array<complex<double>, dim, 1>> ds;

    /**
        The unitary transformation,
        \f$U_0\f$,
        that diagonalises the drift Hamiltonian:
        \f$H_0=U_0D_0U_0^\dagger\f$.

        Also see \ref d0, \ref us_individual, \ref us_inverse_individual, and
        \ref hs.
    */
    Matrix u0;

    /**
        The inverse of the unitary transformation,
        \f$U_0^\dagger\f$,
        that diagonalises the drift Hamiltonian:
        \f$H_0=U_0D_0U_0^\dagger\f$.

        Also see \ref us_individual, \ref d0, and \ref u0_inverse.
    */
    Matrix u0_inverse;

    /**
        The unitary transformations,
        \f$(U_i^\dagger U_{i-1})_{i=1}^{\texttt{length}}\f$,
        from the eigen basis of
        \f$H_{i-1}\f$
        to the eigen basis of
        \f$H_i\f$.

        Also see \ref us_inverse_individual, \ref u0, and \ref d0.
    */
    vector<Matrix> us;

    /**
        The unitary transformations,
        \f$\left(U_i\right)_{i=1}^{\texttt{length}}\f$,
        that diagonalise the control Hamiltonians:
        \f$H_i=U_iD_iU_i^\dagger\f$
        for all
        \f$i\in\left[\texttt{length}\right]\f$.

        Also see \ref us_individual, and \ref us_inverse_individual.
    */
    vector<Matrix> us_individual;

    /**
        The inverse of the unitary transformations,
        \f$(U_i^\dagger)_{i=1}^{\texttt{length}}\f$,
        that diagonalise the control Hamiltonians:
        \f$H_i=U_iD_iU_i^\dagger\f$
        for all
        \f$i\in\left[\texttt{length}\right]\f$.

        Also see \ref us_inverse_individual, \ref us, \ref ds, and \ref hs.
    */
    vector<Matrix> us_inverse_individual;

    /**
        The control Hamiltonians:
        \f$H_i\f$
        for all
        \f$i\in\left[\texttt{length}\right]\f$.

        Also see \ref us_individual, \ref us, \ref ds, and \ref hs.
    */
    vector<Matrix> hs;

    /**
        The unitary transformation,
        \f$U_0^\dagger U_{\texttt{length}}\f$,
        from the eigen basis of
        \f$H_{\texttt{length}}\f$
        to the eigen basis of
        \f$H_0\f$.

        Also see \ref us_individual, \ref us_inverse_individual, and \ref ds.
    */
    Matrix u0_inverse_u_last;
    
    /**
        Initialises a new unitary evolver with the Hamiltonian
        @f[
            H(t)=H_0+\sum_{j=1}^{\texttt{length}}a_j(t)H_j,
        @f]
        where \f$H_0\f$ is the drift Hamiltonian and \f$H_j\f$ are the control
        Hamiltonians modulated by control amplitudes \f$a_j(t)\f$ which need
        not be specified during initialisation.

        @param drift_hamiltonian The drift Hamiltonian.
        @param control_hamiltonians The control Hamiltonians.

        Also see \ref us_individual, \ref us, and \ref u0_inverse.
    */
    template<typename T = Matrix>
    UnitaryEvolver(std::enable_if_t<std::is_same<T, DMatrix<dim, dim>>::value,
                   DMatrix<dim, dim>> drift_hamiltonian,
                   DMatrix<dim_x_n_ctrl, dim> control_hamiltonians
                  ) {
        typedef Eigen::SelfAdjointEigenSolver<DMatrix<dim, dim>> EigSolver;
        EigSolver eigs0(drift_hamiltonian);
        d0 = minus_i*eigs0.eigenvalues().array();
        u0 = eigs0.eigenvectors();
        u0_inverse = u0.adjoint();

        Eigen::Index l = drift_hamiltonian.rows();
        length = control_hamiltonians.rows()/control_hamiltonians.cols();
        
        if (length == 0) {
            u0_inverse_u_last = DMatrix<dim, dim>::Identity(l, l);
            us.push_back(u0);
        } else {
            hs.push_back(
                (DMatrix<dim, dim>)(control_hamiltonians.block(0, 0, l, l)));
            EigSolver eigs(hs.back());
            ds.push_back((Eigen::Array<complex<double>, dim, 1>)(
                minus_i*eigs.eigenvalues().array()));
            DMatrix<dim, dim> u = eigs.eigenvectors();
            us_individual.push_back(u);
            us_inverse_individual.push_back(u.adjoint());
            us.push_back(u.adjoint() * u0);
            us.push_back(u); // append twice
            for (size_t i = 1; i < length; i++) {
                hs.push_back((DMatrix<dim, dim>)(control_hamiltonians.block(
                    i*l, 0, l, l)));
                eigs = EigSolver(hs.back());
                ds.push_back((Eigen::Array<complex<double>, dim, 1>)(
                    minus_i*eigs.eigenvalues().array()));
                u = eigs.eigenvectors();
                us_individual.push_back(u);
                us_inverse_individual.push_back(u.adjoint());
                us.back() = u.adjoint() * us.back();
                us.push_back(u);
            }
            u0_inverse_u_last.noalias() = u0_inverse * u;
        }
    };

    // Sparse initialisation
    /**
        Initialises a new unitary evolver with a sparse Hamiltonian.
        @f[
            H(t)=H_0+\sum_{j=1}^{\texttt{length}}a_j(t)H_j,
        @f]
        where \f$H_0\f$ is the drift Hamiltonian and \f$H_j\f$ are the control
        Hamiltonians modulated by control amplitudes \f$a_j(t)\f$ which need
        not be specified during initialisation.

        @param drift_hamiltonian The drift Hamiltonian.
        @param control_hamiltonians The control Hamiltonians.
    */
    template<typename T = Matrix>
    UnitaryEvolver(std::enable_if_t<std::is_same<T, SMatrix>::value,
                   DMatrix<dim, dim>> drift_hamiltonian,
                   DMatrix<dim_x_n_ctrl, dim> control_hamiltonians) {
        typedef Eigen::SelfAdjointEigenSolver<DMatrix<dim, dim>> EigSolver;
        EigSolver eigs0(drift_hamiltonian);
        d0 = minus_i*eigs0.eigenvalues().array();
        u0 = eigs0.eigenvectors().sparseView();
        u0_inverse = u0.adjoint();

        Eigen::Index l = drift_hamiltonian.rows();
        length = control_hamiltonians.rows()/control_hamiltonians.cols();
        if (length == 0) {
            u0_inverse_u_last = DMatrix<dim, dim>::Identity(l, l).sparseView();
            us.push_back(u0);
        } else {
            hs.push_back(control_hamiltonians.block(0, 0, l, l).sparseView());
            EigSolver eigs((DMatrix<dim, dim>)(hs.back()));
            ds.push_back((Eigen::Array<complex<double>, dim, 1>)(
                                           minus_i*eigs.eigenvalues().array()));
            Matrix u = eigs.eigenvectors().sparseView();
            us_individual.push_back(u);
            us_inverse_individual.push_back(u.adjoint());
            us.push_back((u.adjoint() * u0).pruned());
            us.push_back(u); // append twice
            for (size_t i = 1; i < length; i++) {
                hs.push_back(
                    control_hamiltonians.block(i*l, 0, l, l).sparseView());
                eigs = EigSolver((DMatrix<dim, dim>)(hs.back()));
                ds.push_back((Eigen::Array<complex<double>, dim, 1>)(
                                           minus_i*eigs.eigenvalues().array()));
                u = eigs.eigenvectors().sparseView();
                us_individual.push_back(u);
                us_inverse_individual.push_back(u.adjoint());
                us.back() = (u.adjoint() * us.back()).pruned();
                us.push_back(u);
            }
            u0_inverse_u_last = u0_inverse * u;
        }
    };

    /**
        Initialises a new unitary evolver using the struct attributes.

        @param l The number of control Hamiltonians. Initialises
        \ref length.
        @param d0 The eigenvalues,
        \f$\operatorname{diag}(D_0)\f$,
        of the drift Hamiltonian:
        \f$H_0=U_0D_0U_0^\dagger\f$.
        Initialises \ref d0.
        @param ds The eigenvalues,
        \f$\left(\operatorname{diag}(D_i)\right)_{i=1}^{\texttt{length}}\f$,
        of the control Hamiltonians:
        \f$H_i=U_iD_iU_i^\dagger\f$
        for all
        \f$i\in\left[\texttt{length}\right]\f$.
        Initialises \ref ds.
        @param u0 The unitary transformation,
        \f$U_0\f$,
        that diagonalises the drift Hamiltonian:
        \f$H_0=U_0D_0U_0^\dagger\f$.
        Initialises \ref u0.
        @param u0_inverse The inverse of the unitary transformation,
        \f$U_0^\dagger\f$,
        that diagonalises the drift Hamiltonian:
        \f$H_0=U_0D_0U_0^\dagger\f$.
        Initialises \ref u0_inverse.
        @param us The unitary transformations,
        \f$(U_i^\dagger U_{i-1})_{i=1}^{\texttt{length}}\f$,
        from the eigen basis of
        \f$H_{i-1}\f$
        to the eigen basis of
        \f$H_i\f$.
        Initialises \ref us.
        @param us_individual The unitary transformations,
        \f$\left(U_i\right)_{i=1}^{\texttt{length}}\f$,
        that diagonalise the control Hamiltonians:
        \f$H_i=U_iD_iU_i^\dagger\f$
        for all
        \f$i\in\left[\texttt{length}\right]\f$.
        Initialises \ref us_individual.
        @param us_inverse_individual The inverse of the unitary transformations,
        \f$(U_i^\dagger)_{i=1}^{\texttt{length}}\f$,
        that diagonalise the control Hamiltonians:
        \f$H_i=U_iD_iU_i^\dagger\f$
        for all
        \f$i\in\left[\texttt{length}\right]\f$.
        Initialises \ref us_inverse_individual.
        @param control_hamiltonians The control Hamiltonians:
        \f$H_i\f$
        for all
        \f$i\in\left[\texttt{length}\right]\f$.
        Initialises \ref hs.
        @param u0_inverse_u_last The unitary transformation,
        \f$U_0^\dagger U_{\texttt{length}}\f$,
        from the eigen basis of
        \f$H_{\texttt{length}}\f$
        to the eigen basis of
        \f$H_0\f$.
        Initialises \ref u0_inverse_u_last.
    */
    UnitaryEvolver(size_t l,
                   Eigen::Array<complex<double>, dim, 1> d0,
                   vector<Eigen::Array<complex<double>, dim, 1>> ds,
                   Matrix u0,
                   Matrix u0_inverse,
                   vector<Matrix> us,
                   vector<Matrix> us_individual,
                   vector<Matrix> us_inverse_individual,
                   vector<Matrix> control_hamiltonians,
                   Matrix u0_inverse_u_last)
        : length{l},
          d0{d0},
          ds{ds},
          u0{u0},
          u0_inverse{u0_inverse},
          us{us},
          us_individual{us_individual},
          us_inverse_individual{us_inverse_individual},
          hs{control_hamiltonians},
          u0_inverse_u_last{u0_inverse_u_last}
        {};
        
    /**
        Propagates the state vector using the first-order Suzuki-Trotter
        expansion. More precisely, a state vector,
        \f$\psi(0)\f$,
        is evolved under the differential equation
        @f[
        \dot\psi=-iH\psi
        @f]
        using the first-order Suzuki-Trotter expansion:     
        @f[
        \begin{align}
            \psi(N\Delta t)&=\prod_{i=1}^N\prod_{j=0}^{\texttt{length}}
                e^{-ia_{ij}H_j\Delta t}\psi(0)+\mathcal E\\
            &=\prod_{i=1}^N\prod_{j=0}^{\texttt{length}}
                U_je^{-ia_{ij}D_j\Delta t}U_j^\dagger\psi(0)+\mathcal E.
        \end{align}
        @f]
        where
        \f$a_{nj}\coloneqq a(n\Delta t)\f$,
        we set
        \f$a_{n0}=1\f$
        for notational ease, and the additive error
        \f$\mathcal E\f$
        is
        @f[
        \begin{align}
        \mathcal E&=\mathcal O\left(
            \Delta t^2\left[\sum_{i=1}^N\sum_{j=1}^{\texttt{length}}\dot a_{ij}
            \norm{H_j}
            +\sum_{i=1}^N\sum_{j,k=0}^{\texttt{length}}a_{ij}a_{ik}
            \norm{[H_j,H_k]}\right]
            \right)\\
        &=\mathcal O\left(
            N\Delta t^2\cdot\texttt{length}
            \left[\omega E+\alpha^2E^2\cdot\texttt{length}\right]
            \right)
        \end{align}
        @f]
        where \f$\dot a_{nj}\coloneqq\dot a_j(n\Delta t)\f$ and
        @f[
        \begin{align}
            \omega&\coloneqq\max_{\substack{i\in\left[1,N\right]\\
                j\in\left[1,\texttt{length}\right]}}\left|\dot a_{ij}\right|,\\
            \alpha&\coloneqq\max_{\substack{i\in\left[1,N\right]\\
                j\in\left[0,\texttt{length}\right]}}\left|a_{ij}\right|,\\
            E&\coloneqq\max_{j\in\left[0,\texttt{length}\right]}\norm{H_j}.
        \end{align}
        @f]
        Note the error is quadratic in \f$\Delta t\f$ but linear in \f$N\f$.
        We can also view this as being linear in \f$\Delta t\f$ and linear in
        total evolution time \f$N\Delta t\f$. Additionally, by Nyquist's theorem
        this asymptotic error scaling will not be achieved until the time step
        \f$\Delta t\f$ is smaller than \f$\frac{1}{2\Omega}\f$ where
        \f$\Omega\f$ is the largest energy or frequency in the system.

        @param ctrl_amp \f$\left(a_{ij}\right)\f$ The control amplitudes at each
        time step expressed as an \f$N\times\texttt{length}\f$ matrix where the
        element \f$a_{ij}\f$ corresponds to the control amplitude of the
        \f$j\f$th control Hamiltonian at the \f$i\f$th time step.
        @param state \f$\left[\psi(0)\right]\f$ The state vector to propagate.
        @param dt (\f$\Delta t\f$) The time step to propagate by.
        @return The propagated state vector:
        \f$\psi(N\Delta t)\f$.

        Also see ``propagate_collection()`` and ``propagate_all()``.
    */
    DMatrix<dim, 1> propagate(DMatrix<Dynamic, n_ctrl> ctrl_amp,
                              DMatrix<dim, 1> state,
                              double dt) {
        size_t steps = ctrl_amp.rows();
        DMatrix<dim, 1> exp_d0 = (d0*dt).exp().matrix();

        ctrl_amp *= dt;
        state = u0_inverse * state;
        state = state.cwiseProduct(exp_d0);
        for (size_t j = 0; j < length; j++) {
            state = us[j] * state;
            state = state.cwiseProduct((ds[j]*ctrl_amp(0, j)).exp().matrix());
        }
        for (size_t i = 1; i < steps; i++) {
            state = u0_inverse_u_last * state;
            state = state.cwiseProduct(exp_d0);
            for (size_t j = 0; j < length; j++) {
                state = us[j] * state;
                state = state.cwiseProduct(
                    (ds[j]*ctrl_amp(i, j)).exp().matrix());
            }
        }
        state = us[length] * state;
        return state;
    };

    /**
        Propagates a collection of state vectors using the first-order
        Suzuki-Trotter expansion. More precisely, a collection of state vectors,
        \f$\left(\psi_k(0)\right)_{k}\f$,
        are evolved under the differential equation
        @f[
        \dot\psi_k=-iH\psi_k
        @f]
        using the first-order Suzuki-Trotter expansion:     
        @f[
        \begin{align}
            \psi_k(N\Delta t)&=\prod_{i=1}^N\prod_{j=0}^{\texttt{length}}
                e^{-ia_{ij}H_j\Delta t}\psi_k(0)+\mathcal E\\
            &=\prod_{i=1}^N\prod_{j=0}^{\texttt{length}}
                U_je^{-ia_{ij}D_j\Delta t}U_j^\dagger\psi_k(0)+\mathcal E.
        \end{align}
        @f]
        where
        \f$a_{nj}\coloneqq a(n\Delta t)\f$,
        we set
        \f$a_{n0}=1\f$
        for notational ease, and the addative error
        \f$\mathcal E\f$
        is
        @f[
        \begin{align}
        \mathcal E&=\mathcal O\left(
            \Delta t^2\left[\sum_{i=1}^N\sum_{j=1}^{\texttt{length}}\dot a_{ij}
            \norm{H_j}
            +\sum_{i=1}^N\sum_{j,k=0}^{\texttt{length}}a_{ij}a_{ik}
            \norm{[H_j,H_k]}\right]
            \right)\\
        &=\mathcal O\left(
            N\Delta t^2\cdot\texttt{length}
            \left[\omega E+\alpha^2E^2\cdot\texttt{length}\right]
            \right)
        \end{align}
        @f]
        where \f$\dot a_{nj}\coloneqq\dot a_j(n\Delta t)\f$ and
        @f[
        \begin{align}
            \omega&\coloneqq\max_{\substack{i\in\left[1,N\right]\\
                j\in\left[1,\texttt{length}\right]}}\left|\dot a_{ij}\right|,\\
            \alpha&\coloneqq\max_{\substack{i\in\left[1,N\right]\\
                j\in\left[0,\texttt{length}\right]}}\left|a_{ij}\right|,\\
            E&\coloneqq\max_{j\in\left[0,\texttt{length}\right]}\norm{H_j}.
        \end{align}
        @f]
        Note the error is quadratic in \f$\Delta t\f$ but linear in \f$N\f$.
        We can also view this as being linear in \f$\Delta t\f$ and linear in
        total evolution time \f$N\Delta t\f$. Additionally, by Nyquist's theorem
        this asymptotic error scaling will not be achieved until the time step
        \f$\Delta t\f$ is smaller than \f$\frac{1}{2\Omega}\f$ where
        \f$\Omega\f$ is the largest energy or frequency in the system.

        @tparam l The number of state vectors to propagate.
        @param ctrl_amp \f$\left(a_{ij}\right)\f$ The control amplitudes at each
        time step expressed as an \f$N\times\texttt{length}\f$ matrix where the
        element \f$a_{ij}\f$ corresponds to the control amplitude of the
        \f$j\f$th control Hamiltonian at the \f$i\f$th time step.
        @param states \f$\left[\left(\psi(0)\right)_{k}\right]\f$ A collection
        of state vectors to propagate expressed as a matrix with each column
        corresponding to a state vector.
        @param dt (\f$\Delta t\f$) The time step to propagate by.
        @return The propagated state vectors:
        \f$\left(\psi_k(N\Delta t)\right)_k\f$.

        Also see ``propagate()`` and ``propagate_all()``.
    */
    template<int l = Dynamic>
    DMatrix<dim, l> propagate_collection(DMatrix<Dynamic, n_ctrl> ctrl_amp,
                                         DMatrix<dim, l> states,
                                         double dt) {
        size_t steps = ctrl_amp.rows();
        Eigen::Array<complex<double>, dim, 1> exp_d0 = (d0*dt).exp();

        ctrl_amp *= dt;
        states = u0_inverse*states;
        states = states.array().colwise() * exp_d0;
        for (size_t j = 0; j < length; j++) {
            states = us[j] * states;
            states = states.array().colwise() * (ds[j]*ctrl_amp(0, j)).exp();
        }
        for (size_t i = 1; i < steps; i++) {
            states = u0_inverse_u_last * states;
            states = states.array().colwise() * exp_d0;
            for (size_t j = 0; j < length; j++) {
                states = us[j] * states;
                states = states.array().colwise()
                        *(ds[j]*ctrl_amp(i, j)).exp();
            }
        }
        states = us[length] * states;
        return states;
    };

    /**
        Propagates the state vector using the first-order Suzuki-Trotter
        expansion and returns the resulting state vector at every time step.
        More precisely, a state vector,
        \f$\psi(0)\f$,
        is evolved under the differential equation
        @f[
        \dot\psi=-iH\psi
        @f]
        using the first-order Suzuki-Trotter expansion:     
        @f[
        \begin{align}
            \psi(n\Delta t)&=\prod_{i=1}^n\prod_{j=0}^{\texttt{length}}
                e^{-ia_{ij}H_j\Delta t}\psi(0)+\mathcal E
                \quad\forall n\in\left[0, N\right]\\
            &=\prod_{i=1}^n\prod_{j=0}^{\texttt{length}}
                U_je^{-ia_{ij}D_j\Delta t}U_j^\dagger\psi(0)+\mathcal E
                \quad\forall n\in\left[0, N\right].
        \end{align}
        @f]
        where
        \f$a_{nj}\coloneqq a(n\Delta t)\f$,
        we set
        \f$a_{n0}=1\f$
        for notational ease, and the additive error
        \f$\mathcal E\f$
        is
        @f[
        \begin{align}
        \mathcal E&=\mathcal O\left(
            \Delta t^2\left[\sum_{i=1}^N\sum_{j=1}^{\texttt{length}}\dot a_{ij}
            \norm{H_j}
            +\sum_{i=1}^N\sum_{j,k=0}^{\texttt{length}}a_{ij}a_{ik}
            \norm{[H_j,H_k]}\right]
            \right)\\
        &=\mathcal O\left(
            N\Delta t^2\cdot\texttt{length}
            \left[\omega E+\alpha^2E^2\cdot\texttt{length}\right]
            \right)
        \end{align}
        @f]
        where \f$\dot a_{nj}\coloneqq\dot a_j(n\Delta t)\f$ and
        @f[
        \begin{align}
            \omega&\coloneqq\max_{\substack{i\in\left[1,N\right]\\
                j\in\left[1,\texttt{length}\right]}}\left|\dot a_{ij}\right|,\\
            \alpha&\coloneqq\max_{\substack{i\in\left[1,N\right]\\
                j\in\left[0,\texttt{length}\right]}}\left|a_{ij}\right|,\\
            E&\coloneqq\max_{j\in\left[0,\texttt{length}\right]}\norm{H_j}.
        \end{align}
        @f]
        Note the error is quadratic in \f$\Delta t\f$ but linear in \f$N\f$.
        We can also view this as being linear in \f$\Delta t\f$ and linear in
        total evolution time \f$N\Delta t\f$. Additionally, by Nyquist's theorem
        this asymptotic error scaling will not be achieved until the time step
        \f$\Delta t\f$ is smaller than \f$\frac{1}{2\Omega}\f$ where
        \f$\Omega\f$ is the largest energy or frequency in the system.

        @param ctrl_amp \f$\left(a_{ij}\right)\f$ The control amplitudes at each
        time step expressed as an \f$N\times\texttt{length}\f$ matrix where the
        element \f$a_{ij}\f$ corresponds to the control amplitude of the
        \f$j\f$th control Hamiltonian at the \f$i\f$th time step.
        @param state \f$\left[\psi(0)\right]\f$ The state vector to propagate.
        @param dt (\f$\Delta t\f$) The time step to propagate by.
        @return The propagated state vector at each time step:
        \f$\left(\psi(n\Delta t)\right)_{n=0}^N\f$.

        Also see ``propagate()`` and ``propagate_collection()``.
    */
    DMatrix<dim, Dynamic> propagate_all(DMatrix<Dynamic, n_ctrl> ctrl_amp,
                                        DMatrix<dim, 1> state,
                                        double dt) {
        size_t steps = ctrl_amp.rows();
        DMatrix<dim, 1> exp_d0 = (d0*dt).exp().matrix();
        DMatrix<dim, Dynamic> phi(state.rows(), steps+1);
        ctrl_amp *= dt;
        phi.col(0) = state;
        phi.col(1) = u0_inverse * state;
        phi.col(1) = phi.col(1).cwiseProduct(exp_d0);
        for (size_t j = 0; j < length; j++) {
            phi.col(1) = us[j] * phi.col(1);
            phi.col(1) = phi.col(1).cwiseProduct(
                (ds[j]*ctrl_amp(0, j)).exp().matrix());
        }
        phi.col(1) = us[length] * phi.col(1);
        for (size_t i = 1, k = 2; i < steps; i++, k++) {
            phi.col(k) = u0_inverse * phi.col(i);
            phi.col(k) = phi.col(k).cwiseProduct(exp_d0);
            for (size_t j = 0; j < length; j++) {
                phi.col(k) = us[j] * phi.col(k);
                phi.col(k) = phi.col(k).cwiseProduct(
                    (ds[j]*ctrl_amp(i, j)).exp().matrix());
            }
            phi.col(k) = us[length] * phi.col(k);
        }
        return phi;
    };

    /**
        Calculates the expectation value with respect to an observable of an
        evolved state vector evolved under a control Hamiltonian modulated by
        the control amplitudes. The integration is performed using
        ``propagate()``.

        @param ctrl_amp \f$\left(a_{ij}\right)\f$ The control amplitudes at each
        time step expressed as an \f$N\times\texttt{length}\f$ matrix where the
        element \f$a_{ij}\f$ corresponds to the control amplitude of the
        \f$j\f$th control Hamiltonian at the \f$i\f$th time step.
        @param state \f$\left[\psi(0)\right]\f$ The state vector to propagate.
        @param dt (\f$\Delta t\f$) The time step to propagate by.
        @param observable \f$(\hat O)\f$ The observable to calculate the
        expectation value of.
        @return The expectation value of the observable:
        \f$\langle\hat O\rangle
            \equiv\psi^\dagger(N\Delta t)\hat O\psi(N\Delta t)\f$.

        Also see ``evolved_expectation_value_all()``
    */
    complex<double> evolved_expectation_value(DMatrix<Dynamic, n_ctrl> ctrl_amp,
                                              DMatrix<dim, 1> state,
                                              double dt,
                                              DMatrix<dim, dim> observable) {
        state = propagate(ctrl_amp, state, dt);
        return state.dot(observable * state); // dot conjugates and transposes
    };

    /**
        Calculates the expectation values with respect to an observable of a
        time series of state vectors evolved under a control Hamiltonian
        modulated by the control amplitudes. The integration is performed using
        ``propagate_all()``.

        @param ctrl_amp \f$\left(a_{ij}\right)\f$ The control amplitudes at each
        time step expressed as an \f$N\times\texttt{length}\f$ matrix where the
        element \f$a_{ij}\f$ corresponds to the control amplitude of the
        \f$j\f$th control Hamiltonian at the \f$i\f$th time step.
        @param state \f$\left[\psi(0)\right]\f$ The state vector to propagate.
        @param dt (\f$\Delta t\f$) The time step to propagate by.
        @param observable \f$(\hat O)\f$ The observable to calculate the
        expectation value of.
        @return The expectation value of the observable:
        \f$\left(\psi^\dagger(n\Delta t)\hat O\psi(N\Delta t)\right)_{n=0}^N\f$.

        Also see ``evolved_expectation_value()``
    */
    DMatrix<Dynamic, 1>
    evolved_expectation_value_all(DMatrix<Dynamic, n_ctrl> ctrl_amp,
                                  DMatrix<dim, 1> state,
                                  double dt,
                                  DMatrix<dim, dim> observable) {
        DMatrix<dim, Dynamic> phi = propagate_all(ctrl_amp, state, dt);
        return (phi.adjoint() * observable * phi).diagonal();
    };

    /**
        Calculates the real inner product of an evolved state vector with a
        fixed vector. The evolved state vector is evolved under a control
        Hamiltonian modulated by the control amplitudes. The integration is
        performed using ``propagate()``.

        @param ctrl_amp \f$\left(a_{ij}\right)\f$ The control amplitudes at each
        time step expressed as an \f$N\times\textrm{length}\f$ matrix where the
        element \f$a_{ij}\f$ corresponds to the control amplitude of the
        \f$j\f$th control Hamiltonian at the \f$i\f$th time step.
        @param state \f$\left[\psi(0)\right]\f$ The state vector to propagate.
        @param dt (\f$\Delta t\f$) The time step to propagate by.
        @param fixed_vector \f$(\xi)\f$ The fixed vector to calculate the
        expectation value of.
        @return The inner product of the evolved state vector with the fixed
        vector:
        \f$\sum_{i=1}^\texttt{dim}\xi_i\psi_i(N\Delta t)\f$.

        Also see ``evolved_inner_product_all()``.
    */
    complex<double> evolved_inner_product(DMatrix<Dynamic, n_ctrl> ctrl_amp,
                                          DMatrix<dim, 1> state,
                                          double dt,
                                          DMatrix<1, dim> fixed_vector) {
        state = propagate(ctrl_amp, state, dt);
        return fixed_vector*state;
    };

    /**
        Calculates the real inner products of a time series of evolved state
        vectors with a fixed vector. The evolved state vector is evolved under a
        control Hamiltonian modulated by the control amplitudes. The integration
        is performed using ``propagate_all()``.

        @param ctrl_amp \f$\left(a_{ij}\right)\f$ The control amplitudes at each
        time step expressed as an \f$N\times\textrm{length}\f$ matrix where the
        element \f$a_{ij}\f$ corresponds to the control amplitude of the
        \f$j\f$th control Hamiltonian at the \f$i\f$th time step.
        @param state \f$\left[\psi(0)\right]\f$ The state vector to propagate.
        @param dt (\f$\Delta t\f$) The time step to propagate by.
        @param fixed_vector \f$(\xi)\f$ The fixed vector to calculate the
        expectation value of.
        @return The inner products of the evolved state vectors with the fixed
        vector:
        \f$\left(
        \sum_{i=1}^\texttt{dim}\xi_i\psi_i(n\Delta t)\right)_{n=0}^N\f$.

        Also see ``evolved_inner_product()``.
    */
    DMatrix<Dynamic, 1>
    evolved_inner_product_all(DMatrix<Dynamic, n_ctrl> ctrl_amp,
                              DMatrix<dim, 1> state,
                              double dt,
                              DMatrix<1, dim> fixed_vector) {
        DMatrix<dim, Dynamic> phi = propagate_all(ctrl_amp, state, dt);
        return fixed_vector*phi;
    };

    /**
        Calculates the switching function for a Mayer problem with an
        expectation value as the cost function. More precisely if the cost
        function is
        @f[
        J\left[\vec a(t)\right]\coloneqq\langle\hat O\rangle
            \equiv\psi^\dagger[\vec a(t);T]
            \hat O\psi[\vec a(t);T],
        @f]
        where \f$T=N\Delta t\f$, then the switching function is
        @f[
        \phi_j(t)\coloneqq\frac{\delta J}{\delta a_j(t)}
            =2\operatorname{Im}\left(\psi^\dagger[\vec a(t);T]
            \hat OU(t\to T)H_j\psi[\vec a(t);t]\right).
        @f]
        using the first-order Suzuki-Trotter expansion we can express the
        switching function as
        @f[
        \begin{align}
            &\phi_j(n\Delta t)=\frac{1}{\Delta t}\pdv{J}{a_{nj}}\\
            &=\!2\operatorname{Im}\!\left(\psi^\dagger(T)
                \hat O\!\!\left[\prod_{i>n}^N\prod_{k=1}^{\texttt{length}}
                e^{-ia_{ik}H_k\Delta t}\right]\!\!\!
                \left[\prod_{k=j}^{\texttt{length}}
                e^{-ia_{nk}H_k\Delta t}\right]\!H_j\!\!
                \left[\prod_{k=0}^{j-1}
                e^{-ia_{nk}H_k\Delta t}\right]
                \!\psi(\left[n-1\right]\Delta t)\right),
        \end{align}
        @f]
        where for numerical efficiency we replace
        \f$e^{-ia_{ik}H_k\Delta t}\f$
        with
        \f$U_ke^{-ia_{ik}D_k\Delta t}U_k^\dagger\f$
        as in ``propagate()``.

        @param ctrl_amp \f$\left(a_{ij}\right)\f$ The control amplitudes at each
        time step expressed as an \f$N\times\texttt{length}\f$ matrix where the
        element \f$a_{ij}\f$ corresponds to the control amplitude of the
        \f$j\f$th control Hamiltonian at the \f$i\f$th time step.
        @param state \f$\left[\psi(0)\right]\f$ The initial state vector.
        @param dt (\f$\Delta t\f$) The time step.
        @param cost \f$(\hat O)\f$ The observable to calculate the
        expectation value of.
        @return The expectation value, \f$\psi^\dagger(T)\hat O\psi(T)\f$, and
        the switching function,
        \f$\phi_j(n\Delta t)\f$
        for all
        \f$j\in\left[1,\texttt{length}\right]\f$
        and
        \f$n\in\left[1,N\right]\f$.

        Also see ``gate_switching_function()``.
    */
    std::tuple<complex<double>, Eigen::Matrix<double, Dynamic, n_ctrl>>
    switching_function(DMatrix<Dynamic, n_ctrl> ctrl_amp,
                       DMatrix<dim, 1> state,
                       double dt,
                       DMatrix<dim, dim> cost) {
        // forward propagation
        size_t steps = ctrl_amp.rows();
        DMatrix<dim, 1> exp_d0 = (d0*dt).exp().matrix();
        ctrl_amp *= dt;
        state = u0_inverse * state;
        state = state.cwiseProduct(exp_d0);
        for (size_t j = 0; j < length; j++) {
            state = us[j] * state;
            state = state.cwiseProduct((ds[j]*ctrl_amp(0, j)).exp().matrix());
        }
        for (size_t i = 1; i < steps; i++) {
            state = u0_inverse_u_last * state;
            state = state.cwiseProduct(exp_d0);
            for (size_t j = 0; j < length; j++) {
                state = us[j] * state;
                state = state.cwiseProduct(
                                         (ds[j]*ctrl_amp(i, j)).exp().matrix());
            }
        }
        state = us[length] * state;
        // Lambda
        DMatrix<dim, 1> lambda = cost * state;
        // Energy
        complex<double> E = state.dot(lambda);
        // Back propagation
        DMatrix<dim, 1> exp_d0_c = exp_d0.conjugate();
        ctrl_amp *= -1;
        DMatrix<Dynamic, n_ctrl> out = DMatrix<Dynamic, n_ctrl>::Zero(steps,
                                                                      length);
        for (size_t k = steps-1; true;) {
            for (size_t j = length-1; true;) {
                out(k, j) = lambda.dot(hs[j] * state);
                state = us_inverse_individual[j] * state;
                state = state.cwiseProduct(
                                         (ds[j]*ctrl_amp(k, j)).exp().matrix());
                state = us_individual[j] * state;
                lambda = us_inverse_individual[j] * lambda;
                lambda = lambda.cwiseProduct(
                    (ds[j]*ctrl_amp(k, j)).exp().matrix());
                lambda = us_individual[j] * lambda;
                if (j == 0) {
                    break;
                }
                j--;
            }
            if (k == 0) {
                break;
            }
            state = u0_inverse * state;
            state = state.cwiseProduct(exp_d0_c);
            state = u0 * state;
            lambda = u0_inverse * lambda;
            lambda = lambda.cwiseProduct(exp_d0_c);
            lambda = u0 * lambda;
            k--;
        }
        return
            std::tuple<complex<double>, Eigen::Matrix<double, Dynamic, n_ctrl>>(
                E,
                2*out.imag()
        );
    };
    /**
        Computes the unitary corresponding to the evolution under the
        differential equation
        @f[
        \dot U=-iHU.
        @f]
        The computation is performed using the first-order Suzuki-Trotter
        expansion:     
        @f[
        \begin{align}
            U(N\Delta t)&=\prod_{i=1}^N\prod_{j=0}^{\textrm{length}}
                e^{-ia_{ij}H_j\Delta t}+\mathcal E\\
            &=\prod_{i=1}^N\prod_{j=0}^{\textrm{length}}
                U_je^{-ia_{ij}D_j\Delta t}U_j^\dagger+\mathcal E.
        \end{align}
        @f]
        where
        \f$a_{nj}\coloneqq a(n\Delta t)\f$,
        we set
        $a_{n0}=1$
        for notational ease, and the additive error
        \f$\mathcal E\f$
        is
        @f[
        \begin{align}
        \mathcal E&=\mathcal O\left(
            \Delta t^2\left[\sum_{i=1}^N\sum_{j=1}^{\textrm{length}}\dot a_{ij}
            \norm{H_j}
            +\sum_{i=1}^N\sum_{j,k=0}^{\textrm{length}}a_{ij}a_{ik}
            \norm{[H_j,H_k]}\right]
            \right)\\
        &=\mathcal O\left(
            N\Delta t^2\textrm{length}\left[\omega E+\alpha^2+E^2\right]
            \right)
        \end{align}
        @f]
        where \f$\dot a_{nj}\coloneqq\dot a_j(n\Delta t)\f$ and
        @f[
        \begin{align}
            \omega&\coloneqq\max_{\substack{i\in\left[1,N\right]\\
                j\in\left[1,\textrm{length}\right]}}\left|\dot a_{ij}\right|,\\
            \alpha&\coloneqq\max_{\substack{i\in\left[1,N\right]\\
                j\in\left[0,\textrm{length}\right]}}\left|a_{ij}\right|,\\
            E&\coloneqq\max_{j\in\left[0,\textrm{length}\right]}\norm{H_j}.
        \end{align}
        @f]
        Note the error is quadratic in \f$\Delta t\f$ but linear in \f$N\f$.
        We can also view this as being linear in \f$\Delta t\f$ and linear in
        total evolution time \f$N\Delta t\f$. Additionally, by Nyquist's theorem
        this asymptotic error scaling will not be achieved until the time step
        \f$\Delta t\f$ is smaller than \f$\frac{1}{2\Omega}\f$ where
        \f$\Omega\f$ is the largest energy or frequency in the system.

        @param ctrl_amp \f$\left(a_{ij}\right)\f$ The control amplitudes at each
        time step expressed as an \f$N\times\textrm{length}\f$ matrix where the
        element \f$a_{ij}\f$ corresponds to the control amplitude of the
        \f$j\f$th control Hamiltonian at the \f$i\f$th time step.
        @param dt (\f$\Delta t\f$) The time step to evolve by.
        @return The unitary corresponding to the evolution:
        \f$U(N\Delta t)\f$.

        Also see ``propagate()``, ``propagate_all()``, and
        ``propagate_collection()``.
    */
    DMatrix<dim, dim> get_evolution(DMatrix<Dynamic, n_ctrl> ctrl_amp,
                                    double dt) {
        size_t steps = ctrl_amp.rows();
        Eigen::Array<complex<double>, dim, 1> exp_d0 = (d0*dt).exp();

        ctrl_amp *= dt;
        DMatrix<dim, dim> u = (DMatrix<dim, dim>)u0_inverse;
        u = u.array().colwise() * exp_d0;
        for (size_t j = 0; j < length; j++) {
            u = us[j] * u;
            u = u.array().colwise() * (ds[j]*ctrl_amp(0, j)).exp();
        }
        for (size_t i = 1; i < steps; i++) {
            u = u0_inverse_u_last * u;
            u = u.array().colwise() * exp_d0;
            for (size_t j = 0; j < length; j++) {
                u = us[j] * u;
                u = u.array().colwise() * (ds[j]*ctrl_amp(i, j)).exp();
            }
        }
        u = us[length] * u;
        return u;
    };
    /**
        
        Calculates the gate infidelity with respect to a target gate of the gate
        produced by the control Hamiltonian modulated by the control amplitudes.
        The integration is performed using
        ``get_evolution()``.

        @param ctrl_amp \f$\left(a_{ij}\right)\f$ The control amplitudes at each
        time step expressed as an \f$N\times\textrm{length}\f$ matrix where the
        element \f$a_{ij}\f$ corresponds to the control amplitude of the
        \f$j\f$th control Hamiltonian at the \f$i\f$th time step.
        @param dt (\f$\Delta t\f$) The time step to evolve by.
        @param target The target gate to calculate the infidelity with respect
        to.
        @return The gate infidelity with respect to the target gate:
        @f[
        \mathcal I(U(N\Delta t), \texttt{target})
        \coloneqq 1-\frac{
        \left|\Tr\left[\texttt{target}^\dagger\cdot U(N\Delta t)\right]\right|^2
        +\texttt{dim}}{\texttt{dim}(\texttt{dim}+1)}.
        @f]

        Also see `unitary_gate_infidelity()`.
    */
    double evolved_gate_infidelity(DMatrix<Dynamic, n_ctrl> ctrl_amp,
                                   double dt,
                                   DMatrix<dim, dim> target) {
        DMatrix<dim, dim> U = get_evolution(ctrl_amp, dt);
        return unitary_gate_infidelity(U, target);
    };
    /**
        Calculates the switching function for a Mayer problem with the
        gate infidelity as the cost function. More precisely if the cost
        function is
        @f[
        J\left[\vec a(t)\right]
        \coloneqq\mathcal I(U\left[\vec a(t); T\right], \texttt{target})
        \coloneqq 1-\frac{\left|\Tr\left[\texttt{target}^\dagger
        \cdot U\left[\vec a(t); T\right]\right]\right|^2
        +\texttt{dim}}{\texttt{dim}(\texttt{dim}+1)}.
        @f]
        where \f$T=N\Delta t\f$, then the switching function is
        @f[
        \begin{align}
        &\phi_j(t)\coloneqq\frac{\delta J}{\delta a_j(t)}\\
        &=\frac{2}{\texttt{dim}(\texttt{dim}+1)}\operatorname{Im}\left(
        \Tr\left[U^\dagger(N\Delta t)\cdot\texttt{target}\right]
        \Tr\left[\texttt{target}^\dagger
        \cdot U(t\to T)H_j U[\vec a(t);t]\right]\right).
        \end{align}
        @f]
        Using the first-order Suzuki-Trotter expansion we can express the
        switching function as
        @f[
        \begin{align}
            &\phi_j(n\Delta t)=\frac{1}{\Delta t}\pdv{J}{a_{nj}}\\
            &=\!\frac{2}{\texttt{dim}(\texttt{dim}+1)}\operatorname{Im}
                \!\left(
                \Tr\!\left[U^\dagger(N\Delta t)\cdot\texttt{target}\right]
                \vphantom{[\prod_{k=j}^{\textrm{length}}}\right.\\
                &\left.\cdot\Tr\!\left[\texttt{target}^\dagger\!\cdot\!
                \left[\prod_{i>n}^N\prod_{k=1}^{\textrm{length}}
                e^{-ia_{ik}H_k\Delta t}\right]\!\!\!
                \left[\prod_{k=j}^{\textrm{length}}
                e^{-ia_{nk}H_k\Delta t}\right]\!H_j\!\!
                \left[\prod_{k=0}^{j-1}
                e^{-ia_{nk}H_k\Delta t}\right]
                \! U(\left[n-1\right]\Delta t)\right]\right),
        \end{align}
        @f]
        where for numerical efficiency we replace
        \f$e^{-ia_{ik}H_k\Delta t}\f$
        with
        \f$U_ke^{-ia_{ik}D_k\Delta t}U_k^\dagger\f$
        as in ``get_evolution()``.

        @param ctrl_amp \f$\left(a_{ij}\right)\f$ The control amplitudes at each
        time step expressed as an \f$N\times\textrm{length}\f$ matrix where the
        element \f$a_{ij}\f$ corresponds to the control amplitude of the
        \f$j\f$th control Hamiltonian at the \f$i\f$th time step.
        @param dt (\f$\Delta t\f$) The time step.
        @param cost \f$(\hat O)\f$ The observable to calculate the
        expectation value of.
        @return The gate infidelity,
        \f$I(U\left[\vec a(t); T\right], \texttt{target})\f$ and
        the switching function,
        \f$\phi_j(n\Delta t)\f$
        for all
        \f$j\in\left[1,\textrm{length}\right]\f$
        and
        \f$n\in\left[1,N\right]\f$.

        Also see ``switching_function()``.
    */
    std::tuple<double, Eigen::Matrix<double, Dynamic, n_ctrl>>
    gate_switching_function(DMatrix<Dynamic, n_ctrl> ctrl_amp,
                            double dt,
                            DMatrix<dim, dim> target) {
        // forward propagation
        size_t steps = ctrl_amp.rows();
        Eigen::Array<complex<double>, dim, 1> exp_d0 = (d0*dt).exp();
        ctrl_amp *= dt;
        DMatrix<dim, dim> u = (DMatrix<dim, dim>)u0_inverse;
        u = u.array().colwise() * exp_d0;
        for (size_t j = 0; j < length; j++) {
            u = us[j] * u;
            u = u.array().colwise() * (ds[j]*ctrl_amp(0, j)).exp();
        }
        for (size_t i = 1; i < steps; i++) {
            u = u0_inverse_u_last * u;
            u = u.array().colwise() * exp_d0;
            for (size_t j = 0; j < length; j++) {
                u = us[j] * u;
                u = u.array().colwise() * (ds[j]*ctrl_amp(i, j)).exp();
            }
        }
        u = us[length] * u;
        // Energy
        complex<double> F = std::conj((target.adjoint()*u).trace());
        // Back propagation
        Eigen::Array<complex<double>, dim, 1> exp_d0_c = exp_d0.conjugate();
        ctrl_amp *= -1;
        DMatrix<Dynamic, n_ctrl> out = DMatrix<Dynamic, n_ctrl>::Zero(steps,
                                                                      length);
        for (size_t k = steps-1; true;) {
            for (size_t j = length-1; true;) {
                out(k, j) = (target.adjoint()*hs[j]*u).trace();
                u = us_inverse_individual[j] * u;
                u = u.array().colwise() * (ds[j] * ctrl_amp(k, j)).exp();
                u = us_individual[j] * u;
                target = us_inverse_individual[j] * target;
                target = target.array().colwise()
                         * (ds[j] * ctrl_amp(k, j)).exp();
                target = us_individual[j] * target;
                if (j == 0) {
                    break;
                }
                j--;
            }
            if (k == 0) {
                break;
            }
            u = u0_inverse * u;
            u = u.array().colwise() * exp_d0_c;
            u = u0 * u;
            target = u0_inverse * target;
            target = target.array().colwise() * exp_d0_c;
            target = u0 * target;
            k--;
        }
        Eigen::Index d = target.rows();
        double normalisation = 1.0/(d*(d+1));
        return std::tuple<double, Eigen::Matrix<double, Dynamic, n_ctrl>>(
            1-normalisation*(std::pow(std::abs(F), 2)+d),
            -2*normalisation*(F*out).imag()
        );
    };
};
}