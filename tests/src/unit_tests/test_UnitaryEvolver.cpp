#define _USE_MATH_DEFINES
#include <cmath>
#include <complex>

#include <catch2/catch.hpp>

#include <Suzuki-Trotter-Evolver/UnitaryEvolver.hpp>

using namespace Suzuki_Trotter_Evolver;

using std::complex;

typedef Eigen::Matrix<complex<double>, Dynamic, Dynamic> Matrix;
typedef Eigen::SparseMatrix<complex<double>> SMatrix;
typedef Eigen::Array<complex<double>, Dynamic, Dynamic> Array;

// Initialisation Tests
TEST_CASE("UnitaryEvolver dense Hamiltonian initialisation",
          "[UnitaryEvolver Initialisation]") {
    complex<double> i (0, 1);
    Matrix h0 {{0,  1},  // Pauli X
               {1,  0}}; 
    Matrix hs {{-1,  0}, // Pauli Z 
               { 0, 1},
               { 0, i},  // Pauli Y
               {-i, 0}};

    Matrix H{{-1, 1},
             { 1, 1}};
    H /= std::sqrt(2);
    Matrix SH{{1,  1},
              {i, -i}};
    SH /= std::sqrt(2);
    Array eigs = {{ i},
                  {-i}};

    UnitaryEvolver<Dynamic, Dynamic, Matrix> evolver(h0, hs);

    REQUIRE(evolver.length == 2);

    REQUIRE(evolver.d0.isApprox(eigs));

    REQUIRE(evolver.ds.size() == evolver.length);
    REQUIRE(evolver.ds[0].isApprox(eigs));
    REQUIRE(evolver.ds[1].isApprox(eigs));

    REQUIRE(evolver.u0.isApprox(H)); // X = H Z H where H is the
                                     //  2x2 Hadamard matrix.
    
    REQUIRE(evolver.u0_inverse.isApprox(H));

    REQUIRE(evolver.us.size() == evolver.length + 1);
    REQUIRE(evolver.us[0].isApprox(H));
    REQUIRE(evolver.us[1].isApprox(SH.adjoint()));
    REQUIRE(evolver.us[2].isApprox(SH));

    REQUIRE(evolver.us_individual.size() == evolver.length);
    REQUIRE(evolver.us_individual[0].isApprox(Matrix::Identity(2, 2)));
    REQUIRE(evolver.us_individual[1].isApprox(SH));

    REQUIRE(evolver.us_inverse_individual.size() == evolver.length);
    REQUIRE(evolver.us_inverse_individual[0].isApprox(Matrix::Identity(2, 2)));
    REQUIRE(evolver.us_inverse_individual[1].isApprox(SH.adjoint()));

    REQUIRE(evolver.hs.size() == evolver.length);
    REQUIRE(evolver.hs[0].isApprox(hs.block(0, 0, 2, 2)));
    REQUIRE(evolver.hs[1].isApprox(hs.block(2, 0, 2, 2)));

    REQUIRE(evolver.u0_inverse_u_last.isApprox(H * SH));
}

TEST_CASE("UnitaryEvolver dense Hamiltonian initialisation no controls",
          "[UnitaryEvolver Initialisation]") {
complex<double> i (0, 1);
Matrix h0 {{0,  1},  // Pauli X
           {1,  0}}; 
Matrix hs = Matrix::Zero(0, 2);

Matrix H{{-1, 1},
         { 1, 1}};
H /= std::sqrt(2);
Array eigs = {{ i},
              {-i}};

UnitaryEvolver<Dynamic, Dynamic, Matrix> evolver(h0, hs);

REQUIRE(evolver.length == 0);

REQUIRE(evolver.d0.isApprox(eigs));

REQUIRE(evolver.ds.size() == evolver.length);

REQUIRE(evolver.u0.isApprox(H)); // X = H Z H where H is the
                                 //  2x2 Hadamard matrix.

REQUIRE(evolver.u0_inverse.isApprox(H));

REQUIRE(evolver.us.size() == evolver.length + 1);
REQUIRE(evolver.us[0].isApprox(H));

REQUIRE(evolver.us_individual.size() == evolver.length);

REQUIRE(evolver.us_inverse_individual.size() == evolver.length);

REQUIRE(evolver.hs.size() == evolver.length);

REQUIRE(evolver.u0_inverse_u_last.isApprox(Matrix::Identity(2, 2)));
}


TEST_CASE("UnitaryEvolver sparse Hamiltonian initialisation",
          "[UnitaryEvolver Initialisation]") {
complex<double> i (0, 1);
Matrix h0d {{0,  1},  // Pauli X
            {1,  0}}; 
Matrix hsd {{-1,  0}, // Pauli Z 
            { 0, 1},
            { 0, i},  // Pauli Y
            {-i, 0}};
SMatrix h0 = h0d.sparseView();
SMatrix hs = hsd.sparseView();

Matrix H{{-1, 1},
         { 1, 1}};
H /= std::sqrt(2);
Matrix SH{{1,  1},
          {i, -i}};
SH /= std::sqrt(2);
Array eigs = {{ i},
              {-i}};

UnitaryEvolver<Dynamic, Dynamic, SMatrix> evolver(h0, hs);

REQUIRE(evolver.length == 2);

REQUIRE(evolver.d0.isApprox(eigs));

REQUIRE(evolver.ds.size() == evolver.length);
REQUIRE(evolver.ds[0].isApprox(eigs));
REQUIRE(evolver.ds[1].isApprox(eigs));

REQUIRE(evolver.u0.isApprox(H)); // X = H Z H where H is the
                               //  2x2 Hadamard matrix.

REQUIRE(evolver.u0_inverse.isApprox(H));

REQUIRE(evolver.us.size() == evolver.length + 1);
REQUIRE(evolver.us[0].isApprox(H));
REQUIRE(evolver.us[1].isApprox(SH.adjoint()));
REQUIRE(evolver.us[2].isApprox(SH));

REQUIRE(evolver.us_individual.size() == evolver.length);
REQUIRE(evolver.us_individual[0].isApprox(Matrix::Identity(2, 2)));
REQUIRE(evolver.us_individual[1].isApprox(SH));

REQUIRE(evolver.us_inverse_individual.size() == evolver.length);
REQUIRE(evolver.us_inverse_individual[0].isApprox(Matrix::Identity(2, 2)));
REQUIRE(evolver.us_inverse_individual[1].isApprox(SH.adjoint()));

REQUIRE(evolver.hs.size() == evolver.length);
REQUIRE(evolver.hs[0].isApprox(hs.block(0, 0, 2, 2)));
REQUIRE(evolver.hs[1].isApprox(hs.block(2, 0, 2, 2)));

REQUIRE(evolver.u0_inverse_u_last.isApprox(H * SH));
}

TEST_CASE("UnitaryEvolver sparse Hamiltonian initialisation no controls",
          "[UnitaryEvolver Initialisation]") {
complex<double> i (0, 1);
Matrix h0d {{0,  1},  // Pauli X
           {1,  0}}; 
Matrix hsd = Matrix::Zero(0, 2);

SMatrix h0 = h0d.sparseView();
SMatrix hs = hsd.sparseView();

Matrix H{{-1, 1},
         { 1, 1}};
H /= std::sqrt(2);
Array eigs = {{ i},
              {-i}};

UnitaryEvolver<Dynamic, Dynamic, SMatrix> evolver(h0, hs);

REQUIRE(evolver.length == 0);

REQUIRE(evolver.d0.isApprox(eigs));

REQUIRE(evolver.ds.size() == evolver.length);

REQUIRE(evolver.u0.isApprox(H)); // X = H Z H where H is the
                                 //  2x2 Hadamard matrix.

REQUIRE(evolver.u0_inverse.isApprox(H));

REQUIRE(evolver.us.size() == evolver.length + 1);
REQUIRE(evolver.us[0].isApprox(H));

REQUIRE(evolver.us_individual.size() == evolver.length);

REQUIRE(evolver.us_inverse_individual.size() == evolver.length);

REQUIRE(evolver.hs.size() == evolver.length);

REQUIRE(evolver.u0_inverse_u_last.isApprox(Matrix::Identity(2, 2)));
}

TEST_CASE("UnitaryEvolver property initialisation",
          "[UnitaryEvolver Initialisation]") {
complex<double> i (0, 1);
Matrix h0 {{0,  1},  // Pauli X
           {1,  0}}; 
Matrix hs {{-1,  0}, // Pauli Z 
           { 0, 1},
           { 0, i},  // Pauli Y
           {-i, 0}};

Matrix H{{-1, 1},
         { 1, 1}};
H /= std::sqrt(2);
Matrix SH{{1,  1},
          {i, -i}};
SH /= std::sqrt(2);
Array eigs = {{ i},
              {-i}};

UnitaryEvolver<Dynamic, Dynamic, Matrix> ev(h0, hs);
UnitaryEvolver<Dynamic, Dynamic, Matrix> evolver(ev.length,
                                                 ev.d0,
                                                 ev.ds,
                                                 ev.u0,
                                                 ev.u0_inverse,
                                                 ev.us,
                                                 ev.us_individual,
                                                 ev.us_inverse_individual,
                                                 ev.hs,
                                                 ev.u0_inverse_u_last);

REQUIRE(evolver.length == 2);

REQUIRE(evolver.d0.isApprox(eigs));

REQUIRE(evolver.ds.size() == evolver.length);
REQUIRE(evolver.ds[0].isApprox(eigs));
REQUIRE(evolver.ds[1].isApprox(eigs));

REQUIRE(evolver.u0.isApprox(H)); // X = H Z H where H is the
                               //  2x2 Hadamard matrix.

REQUIRE(evolver.u0_inverse.isApprox(H));

REQUIRE(evolver.us.size() == evolver.length + 1);
REQUIRE(evolver.us[0].isApprox(H));
REQUIRE(evolver.us[1].isApprox(SH.adjoint()));
REQUIRE(evolver.us[2].isApprox(SH));

REQUIRE(evolver.us_individual.size() == evolver.length);
REQUIRE(evolver.us_individual[0].isApprox(Matrix::Identity(2, 2)));
REQUIRE(evolver.us_individual[1].isApprox(SH));

REQUIRE(evolver.us_inverse_individual.size() == evolver.length);
REQUIRE(evolver.us_inverse_individual[0].isApprox(Matrix::Identity(2, 2)));
REQUIRE(evolver.us_inverse_individual[1].isApprox(SH.adjoint()));

REQUIRE(evolver.hs.size() == evolver.length);
REQUIRE(evolver.hs[0].isApprox(hs.block(0, 0, 2, 2)));
REQUIRE(evolver.hs[1].isApprox(hs.block(2, 0, 2, 2)));

REQUIRE(evolver.u0_inverse_u_last.isApprox(H * SH));
}

TEST_CASE("UnitaryEvolver property initialisation no controls",
          "[UnitaryEvolver Initialisation]") {
complex<double> i (0, 1);
Matrix h0 {{0,  1},  // Pauli X
           {1,  0}}; 
Matrix hs = Matrix::Zero(0, 2);

Matrix H{{-1, 1},
         { 1, 1}};
H /= std::sqrt(2);
Array eigs = {{ i},
              {-i}};

UnitaryEvolver<Dynamic, Dynamic, Matrix> ev(h0, hs);
UnitaryEvolver<Dynamic, Dynamic, Matrix> evolver(ev.length,
                                                 ev.d0,
                                                 ev.ds,
                                                 ev.u0,
                                                 ev.u0_inverse,
                                                 ev.us,
                                                 ev.us_individual,
                                                 ev.us_inverse_individual,
                                                 ev.hs,
                                                 ev.u0_inverse_u_last);

REQUIRE(evolver.length == 0);

REQUIRE(evolver.d0.isApprox(eigs));

REQUIRE(evolver.ds.size() == evolver.length);

REQUIRE(evolver.u0.isApprox(H)); // X = H Z H where H is the
                                 //  2x2 Hadamard matrix.

REQUIRE(evolver.u0_inverse.isApprox(H));

REQUIRE(evolver.us.size() == evolver.length + 1);
REQUIRE(evolver.us[0].isApprox(H));

REQUIRE(evolver.us_individual.size() == evolver.length);

REQUIRE(evolver.us_inverse_individual.size() == evolver.length);

REQUIRE(evolver.hs.size() == evolver.length);

REQUIRE(evolver.u0_inverse_u_last.isApprox(Matrix::Identity(2, 2)));
}



// Evolution Tests
TEST_CASE("propagate identity", "[UnitaryEvolver Propagation]") {
    complex<double> i (0, 1);
    Matrix h0 = Matrix::Zero(2, 2);
    Matrix hs {{-1,  0}, // Pauli Z 
               { 0, 1},
               { 0, i},  // Pauli Y
               {-i, 0}};

    UnitaryEvolver<Dynamic, Dynamic, Matrix> evolver(h0, hs);

    Matrix ctrl_amp = Matrix::Zero(100, 2);
    double dt = 0.1;

    Matrix initial_state {{1},
                          {0}};

    Matrix output_state = evolver.propagate(ctrl_amp, initial_state, dt);

    REQUIRE(output_state.isApprox(initial_state));
}

TEST_CASE("propagate constant evolution", "[UnitaryEvolver Propagation]") {
    complex<double> i (0, 1);
    Matrix h0 {{1,  0},
               {0, -1}}; 
    Matrix hs {{-1, 0},  // Pauli Z 
               { 0, 1},
               { 0, i},  // Pauli Y
               {-i, 0}};

    UnitaryEvolver<Dynamic, Dynamic, Matrix> evolver(h0, hs);

    Matrix ctrl_amp = Matrix::Zero(100, 2);
    double dt = 0.1;

    Matrix initial_state {{1},
                          {0}};

    Matrix output_state = evolver.propagate(ctrl_amp, initial_state, dt);

    REQUIRE(output_state.isApprox(Matrix{{std::exp(-i*dt*100.)},
                                         {0}}));
}

TEST_CASE("propagate no control Hamiltonians", "[UnitaryEvolver Propagation]") {
    complex<double> i (0, 1);
    Matrix h0 {{1,  0},
               {0,  -1}};
    Matrix hs = Matrix::Zero(0, 2);

    UnitaryEvolver<Dynamic, Dynamic, Matrix> evolver(h0, hs);

    Matrix ctrl_amp = Matrix::Zero(100, 0);
    double dt = 0.1;

    Matrix initial_state {{1},
                          {0}};

    Matrix output_state = evolver.propagate(ctrl_amp, initial_state, dt);

    REQUIRE(output_state.isApprox(Matrix{{std::exp(-i*dt*100.)},
                                         {0}}));
}

TEST_CASE("propagate Rabi oscillation", "[UnitaryEvolver Propagation]") {
    // This test takes the rotating wave approximation in the analytics and so
    //  isApprox is vital.
    complex<double> i (0, 1);
    Matrix h0 {{1,  0},   // Pauli Z 
               {0, -1}}; 
    Matrix hs {{0, 1},    // Pauli X
               {1, 0}};

    UnitaryEvolver<Dynamic, Dynamic, Matrix> evolver(h0, hs);

    
    double dt = 0.001;
    double omega = 2*M_PI*0.001;

    Matrix initial_state {{1},
                          {0}};
    Matrix ctrl_amp = Matrix::Zero(1000000, 1);
    for (size_t k = 0; k < 1000000; k++) {
        ctrl_amp(k, 0) = omega*std::cos(2*k*dt); // Resonant drive
    }
    for (size_t j = 100000; j <= 1000000; j+=100000) {
        Matrix output_state =
            evolver.propagate(ctrl_amp(Eigen::seqN(0,j), Eigen::all),
                              initial_state,
                              dt);
        Matrix rwa_analytic_output_state
            {{   std::exp(-i*(double)j*dt)*std::cos(omega*(double)j*dt/2)},
             {-i*std::exp( i*(double)j*dt)*std::sin(omega*(double)j*dt/2)}};
        REQUIRE(output_state.isApprox(rwa_analytic_output_state, 5e-3));
    }
}


TEST_CASE("propagate_collection identity", "[UnitaryEvolver Propagation]") {
    complex<double> i (0, 1);
    Matrix h0 = Matrix::Zero(2, 2);
    Matrix hs {{-1,  0}, // Pauli Z 
               { 0, 1},
               { 0, i},  // Pauli Y
               {-i, 0}};

    UnitaryEvolver<Dynamic, Dynamic, Matrix> evolver(h0, hs);

    Matrix ctrl_amp = Matrix::Zero(100, 2);
    double dt = 0.1;

    Matrix initial_states {{1, 0},
                           {0, 1}};

    Matrix output_state = evolver.propagate_collection(ctrl_amp,
                                                       initial_states,
                                                       dt);

    REQUIRE(output_state.isApprox(initial_states));
}

TEST_CASE("propagate_collection constant evolution",
          "[UnitaryEvolver Propagation]") {
    complex<double> i (0, 1);
    Matrix h0 {{1,  0},
               {0,  -1}}; 
    Matrix hs {{-1,  0}, // Pauli Z 
               { 0, 1},
               { 0, i},  // Pauli Y
               {-i, 0}};

    UnitaryEvolver<Dynamic, Dynamic, Matrix> evolver(h0, hs);

    Matrix ctrl_amp = Matrix::Zero(100, 2);
    double dt = 0.1;

    Matrix initial_states {{1, 0},
                           {0, 1}};

    Matrix output_state = evolver.propagate_collection(ctrl_amp,
                                                       initial_states,
                                                       dt);

    REQUIRE(output_state.isApprox(Matrix{{std::exp(-i*dt*100.), 0},
                                         {0, std::exp( i*dt*100.)}}));
}

TEST_CASE("propagate_collection no control Hamiltonians",
          "[UnitaryEvolver Propagation]") {
    complex<double> i (0, 1);
    Matrix h0 {{1,  0},
               {0,  -1}};
    Matrix hs = Matrix::Zero(0, 2);

    UnitaryEvolver<Dynamic, Dynamic, Matrix> evolver(h0, hs);

    Matrix ctrl_amp = Matrix::Zero(100, 0);
    double dt = 0.1;

    Matrix initial_states {{1, 0},
                           {0, 1}};

    Matrix output_state = evolver.propagate_collection(ctrl_amp,
                                                       initial_states,
                                                       dt);

    REQUIRE(output_state.isApprox(Matrix{{std::exp(-i*dt*100.), 0},
                                         {0, std::exp( i*dt*100.)}}));
}

TEST_CASE("propagate_collection Rabi oscillation",
          "[UnitaryEvolver Propagation]") {
    Matrix h0 {{1,  0},   // Pauli Z 
               {0, -1}}; 
    Matrix hs {{0, 1},    // Pauli X
               {1, 0}};

    UnitaryEvolver<Dynamic, Dynamic, Matrix> evolver(h0, hs);

    
    double dt = 0.001;
    double omega = 2*M_PI*0.1;


    Matrix initial_states {{1, 0},
                           {0, 1}};
    Matrix ctrl_amp = Matrix::Zero(10000, 1);
    for (size_t k = 0; k < 10000; k++) {
        ctrl_amp(k, 0) = omega*std::cos(2*k*dt); // Resonant drive
    }
    for (size_t j = 1000; j <= 10000; j+=1000) {
        Matrix cropped_ctrl_amp = ctrl_amp(Eigen::seqN(0,j), Eigen::all);
        Matrix output_state = evolver.propagate_collection(cropped_ctrl_amp,
                                                           initial_states,
                                                           dt);
        Matrix propagate_output(2, 2);
        propagate_output << evolver.propagate(cropped_ctrl_amp,
                                              initial_states.col(0),
                                              dt),
                            evolver.propagate(cropped_ctrl_amp,
                                              initial_states.col(1),
                                              dt);
        REQUIRE(output_state.isApprox(propagate_output));
    }
}


TEST_CASE("propagate_all identity", "[UnitaryEvolver Propagation]") {
    complex<double> i (0, 1);
    Matrix h0 = Matrix::Zero(2, 2);
    Matrix hs {{-1,  0}, // Pauli Z 
               { 0, 1},
               { 0, i},  // Pauli Y
               {-i, 0}};

    UnitaryEvolver<Dynamic, Dynamic, Matrix> evolver(h0, hs);

    Matrix ctrl_amp = Matrix::Zero(100, 2);
    double dt = 0.1;

    Matrix initial_state {{1},
                          {0}};

    Matrix output_states = evolver.propagate_all(ctrl_amp, initial_state, dt);

    Matrix target_output_states(2, 101);
    target_output_states << Matrix::Ones(1, 101),
                            Matrix::Zero(1, 101);
    REQUIRE(output_states.isApprox(target_output_states));
}

TEST_CASE("propagate_all constant evolution", "[UnitaryEvolver Propagation]") {
    complex<double> i (0, 1);
    Matrix h0 {{1,  0},
               {0,  -1}}; 
    Matrix hs {{-1,  0}, // Pauli Z 
               { 0, 1},
               { 0, i},  // Pauli Y
               {-i, 0}};

    UnitaryEvolver<Dynamic, Dynamic, Matrix> evolver(h0, hs);

    Matrix ctrl_amp = Matrix::Zero(100, 2);
    double dt = 0.1;

    Matrix initial_state {{1},
                          {0}};

    Matrix output_states = evolver.propagate_all(ctrl_amp, initial_state, dt);

    Matrix target_output_states = Matrix::Zero(2, 101);
    for (size_t k = 0; k < 101; k++) {
        target_output_states(0, k) = std::exp(-i*dt*(double)k);
    }

    REQUIRE(output_states.isApprox(target_output_states));
}

TEST_CASE("propagate_all no control Hamiltonians",
          "[UnitaryEvolver Propagation]") {
    complex<double> i (0, 1);
    Matrix h0 {{1,  0},
               {0,  -1}};
    Matrix hs = Matrix::Zero(0, 2);

    UnitaryEvolver<Dynamic, Dynamic, Matrix> evolver(h0, hs);

    Matrix ctrl_amp = Matrix::Zero(100, 0);
    double dt = 0.1;

    Matrix initial_state {{1},
                          {0}};

    Matrix output_states = evolver.propagate_all(ctrl_amp, initial_state, dt);

    Matrix target_output_states = Matrix::Zero(2, 101);
    for (size_t k = 0; k < 101; k++) {
        target_output_states(0, k) = std::exp(-i*dt*(double)k);
    }

    REQUIRE(output_states.isApprox(target_output_states));
}

TEST_CASE("propagate_all Rabi oscillation",
          "[UnitaryEvolver Propagation]") {
    Matrix h0 {{1,  0},   // Pauli Z 
               {0, -1}}; 
    Matrix hs {{0, 1},    // Pauli X
               {1, 0}};

    UnitaryEvolver<Dynamic, Dynamic, Matrix> evolver(h0, hs);

    
    double dt = 0.01;
    double omega = 2*M_PI;

    Matrix ctrl_amp = Matrix::Zero(100, 1);
    for (size_t k = 0; k < 100; k++) {
        ctrl_amp(k, 0) = omega*std::cos(2*k*dt); // Resonant drive
    }

    Matrix initial_state {{1},
                          {0}};

    Matrix output_states = evolver.propagate_all(ctrl_amp, initial_state, dt);
    Matrix target_output_states(2, 101);
    target_output_states.col(0) = initial_state;

    for (size_t j = 1; j <= 100; j++) {
        target_output_states.col(j) =
            evolver.propagate(ctrl_amp(Eigen::seqN(0,j), Eigen::all),
                              initial_state,
                              dt);
    }
    REQUIRE(output_states.isApprox(target_output_states));
}


// Test Expectation Values
TEST_CASE("evolved_expectation_value", "[UnitaryEvolver Expectation_Value]") {
    complex<double> i (0, 1);
    Matrix h0 {{1,  0},   // Pauli Z 
               {0, -1}}; 
    Matrix hs {{0, 1},    // Pauli X
               {1, 0}};
    Matrix Y  {{0, -i},
               {i,  0}}; // Pauli Y

    UnitaryEvolver<Dynamic, Dynamic, Matrix> evolver(h0, hs);

    
    double dt = 0.001;
    double omega = 2*M_PI*0.1;

    Matrix initial_state {{1},
                          {0}};
    Matrix ctrl_amp = Matrix::Zero(2500, 1);
    for (size_t k = 0; k < 2500; k++) {
        ctrl_amp(k, 0) = omega*std::cos(2*k*dt); // Resonant drive
    }
    complex<double> expval = evolver.evolved_expectation_value(ctrl_amp,
                                                               initial_state,
                                                               dt,
                                                               Y);
    
    Matrix state = evolver.propagate(ctrl_amp, initial_state, dt);
    complex<double> target_expval = (state.adjoint() * Y * state)(0, 0);

    REQUIRE(std::abs(expval - target_expval) <= 1e-15);
}

TEST_CASE("evolved_expectation_value_all",
          "[UnitaryEvolver Expectation_Value]") {
    complex<double> i (0, 1);
    Matrix h0 {{1,  0},   // Pauli Z 
               {0, -1}}; 
    Matrix hs {{0, 1},    // Pauli X
               {1, 0}};
    Matrix Y  {{0, -i},
               {i,  0}}; // Pauli Y

    UnitaryEvolver<Dynamic, Dynamic, Matrix> evolver(h0, hs);

    
    double dt = 0.001;
    double omega = 2*M_PI*0.1;

    Matrix initial_state {{1},
                          {0}};
    Matrix ctrl_amp = Matrix::Zero(2500, 1);
    for (size_t k = 0; k < 2500; k++) {
        ctrl_amp(k, 0) = omega*std::cos(2*k*dt); // Resonant drive
    }
    Matrix expval = evolver.evolved_expectation_value_all(ctrl_amp,
                                                          initial_state,
                                                          dt,
                                                          Y);
    
    Matrix state = evolver.propagate_all(ctrl_amp, initial_state, dt);
    Matrix target_expvals(2501, 1);
    for (size_t k = 0; k < 2501; k++) {
        target_expvals(k, 0) = (state.col(k).adjoint()*Y*state.col(k))(0, 0);
    }

    REQUIRE(expval.isApprox(target_expvals));
}

// Test Switching Function
TEST_CASE("switching_function", "[UnitaryEvolver Switching_Function]") {
    complex<double> i (0, 1);
    Matrix h0 {{1,  0},   // Pauli Z 
               {0, -1}}; 
    Matrix hs {{0, 1},    // Pauli X
               {1, 0},
               {0, -i},   // Pauli Y
               {i,  0}};
    Matrix H {{1,  1},
              {1, -1}};

    UnitaryEvolver<Dynamic, Dynamic, Matrix> evolver(h0, hs);

    double dt = 0.1;
    double omega = 2*M_PI*0.1;

    Matrix initial_state {{1},
                          {0}};
    Matrix ctrl_amp = Matrix::Zero(100, 2);
    for (size_t j = 0; j < 100; j++) {
        ctrl_amp(j, 0) = omega*std::cos(2*j*dt); // Resonant drive
    }

    std::tuple<complex<double>, Matrix> output =
        evolver.switching_function(ctrl_amp, initial_state, dt, H);

    complex<double> expval1 = evolver.evolved_expectation_value(ctrl_amp,
                                                                initial_state,
                                                                dt,
                                                                H);
    REQUIRE(std::abs(std::get<0>(output) - expval1) <= 1e-15);
    double eps = 1e-6;
    Matrix fd_expval(100, 2);
    for (size_t j = 0; j < 100; j++) {
        for (size_t k = 0; k < 2; k++) {
            Matrix new_ctrl_amp = ctrl_amp;
            new_ctrl_amp(j, k) += eps;
            complex<double> expval2 =
                evolver.evolved_expectation_value(new_ctrl_amp,
                                                initial_state,
                                                dt,
                                                H);
            fd_expval(j, k) = (expval2 - expval1) / eps;
        }
    }
    fd_expval /= dt;
    REQUIRE(std::get<1>(output).isApprox(fd_expval, eps));
}