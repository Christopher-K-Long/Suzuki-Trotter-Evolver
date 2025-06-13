// examples/Rabi_oscillation.cpp
#include <cmath>
#include <iostream>

#include <Eigen/Dense>

#include <Suzuki-Trotter-Evolver/UnitaryEvolver.hpp>

using namespace Suzuki_Trotter_Evolver;

using Eigen::MatrixXcd;

// Defining the Pauli matrices
const MatrixXcd X {{0,  1},  // Pauli X
                   {1,  0}};
const MatrixXcd Z {{1,  0},   // Pauli Z
                   {0, -1}};

int main() {
    // Parameters
    double v = 1;
    double f = 1;
    double T = 1;
    int N = 1000;
    MatrixXcd initial_state {{1},
                             {0}};

    // Setting up the evolver
    MatrixXcd drift_hamiltonian = v/2 * Z;
    MatrixXcd control_hamiltonian = X;

    UnitaryEvolver evolver(drift_hamiltonian, control_hamiltonian);

    // Specifying the control amplitudes
    double dt = T / N;
    MatrixXcd ctrl_amp(N, 1);
    
    for (int n = 0; n < N; n++) {
        ctrl_amp(n, 0) = f*std::cos(v * n * dt);
    }

    // Propagating the initial state
    MatrixXcd psi = evolver.propagate(ctrl_amp, initial_state, dt);

    // Outputting the final state
    std::cout<<psi<<std::endl;

    return 0;
}