// examples/compute_switching_function.cpp
#include <cmath>
#include <iostream>

#include <Eigen/Dense>

#include <Suzuki-Trotter-Evolver/UnitaryEvolver.hpp>

using namespace Suzuki_Trotter_Evolver;

using Eigen::MatrixXcd;

// Defining the Pauli matrices
const complex<double> i(0, 1);
const MatrixXcd X {{0,  1},  // Pauli X
                   {1,  0}};
const MatrixXcd Y {{0, -i},  // Pauli Y
                   {i,  0}};
const MatrixXcd Z {{1,  0},   // Pauli Z
                   {0, -1}};

int main() {
    // Parameters
    double v = 1;
    double f = 1;
    double T = 1;
    int N = 20;
    MatrixXcd initial_state {{1},
                             {0}};
    MatrixXcd O {{1,  1},
                 {1, -1}};

    // Setting up the evolver
    MatrixXcd drift_hamiltonian = v/2 * Z;
    MatrixXcd control_hamiltonians(4, 2);
    control_hamiltonians << X, Y;

    UnitaryEvolver evolver(drift_hamiltonian, control_hamiltonians);

    // Specifying the control amplitudes
    double dt = T / N;
    MatrixXcd ctrl_amp(N, 2);
    
    for (int n = 0; n < N; n++) {
        ctrl_amp(n, 0) = f*std::cos(v * n * dt);
    }

    // Computing the expectation value and the switching function
    std::tuple<complex<double>, MatrixXcd> output =
        evolver.switching_function(ctrl_amp, initial_state, dt, O);

    // Unpacking results
    complex<double> expectation_value = std::get<0>(output);
    MatrixXcd switching_function = std::get<1>(output);

    // Outputting the results
    std::cout<<"Expectation value: "<<expectation_value.real()<<std::endl;
    std::cout<<"Switching function:"<<std::endl<<switching_function<<std::endl;

    return 0;
}