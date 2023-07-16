/**
 * @code quantum_dynamics_1d.c
 * @brief calc 1 Dimension Quantum Dynamics with Crank-Nicolson method
 */

#include "quantum_dynamics_1d.h"
#include "SparseCore"
#include "SparseQR"
#include <complex>
#include <QVector>

/**
 * @brief Kronecker delta
 *
 * @param[in] i arg1
 * @param[in] j arg2
 *
 * @return 1:i=j, 0:i!=j
 */
static int delta(const int &i, const int &j)
{
    return ( i == j ) ? 1 : 0;
}

/**
 * @brief potential well
 *
 * @param[in] x
 * @param[in] witdh width of well.
 * @param[in] strength strength of well. strength > 0 means it is hight. strength < 0 means it is depth.
 * @param[in] otherParameter NOT use here.
 *
 * @return potential at x.
 */
std::complex<double> potential_well(const double &x, const double &width,
                                    const std::complex<double> &strength, const QVector<std::complex<double>> &otherParameter)
{
    return ( abs(x) < width ) ? strength : 0;
}

/**
 * @brief woods_saxon potential
 *
 * @param[in] x
 * @param[in] width nuclear radius.
 * @param[in] strength strength of potential. strength > 0 means it is hight. strength < 0 means it is depth.
 * @param[in] otherParameter [0]:a.
 *
 * @return optential at x.
 */
std::complex<double> woods_saxon(const double &x, const double &width,
                                 const std::complex<double> &strength, const QVector<std::complex<double>> &otherParameter)
{
    return strength / (1.0 + exp((abs(x) - width)/otherParameter[0]));
}

/**
 * @brief 2 Range Gausian potential
 *
 * @param[in] x
 * @param[in] width NOT use here.
 * @param[in] strength V0
 * @param[in] otherParamter [1]:V1, [2]:mu0, [3]:mu1

 * @return potential at x.
 */
std::complex<double> two_range_gausian(const double &x, const double &width,
                                       const std::complex<double> &strength, const QVector<std::complex<double>> &otherParameter)
{
    return strength * exp( -otherParameter[2] * x * x ) + otherParameter[1] * exp( -otherParameter[3] * x * x );
}

/**
 * @brief set Solver for 1 Direction Quantum Dynamics
 *
 * @param[out] solver_prev solver for dt < 0
 * @param[out] solver_next solver for dt > 0
 * @param[out] Hamiltonian_1 coefficient matrix for psi(t)
 * @param[out] Hamiltonian_2 coefficient matrix for psi(t + dt)
 * @param[in] div_x stride of x
 * @param[in] x_min min of xAxis
 * @param[in] x_max max of xAxis
 * @param[in] dt stride of time [s・c = fm]
 * @param[in] mass
 * @param[in] width potential parameter
 * @param[in] stength potential parameter
 * @param[in] otherParameter potential parameter
 * @param[in] potential function pointer for potential
 *
 * @return true:success, false:failure
 */
bool set_Solver(Eigen::SparseQR< Eigen::SparseMatrix<Eigen::dcomplex>, Eigen::COLAMDOrdering<int> > &solver_prev,
               Eigen::SparseQR< Eigen::SparseMatrix<Eigen::dcomplex>, Eigen::COLAMDOrdering<int> > &solver_next,
               Eigen::SparseMatrix< Eigen::dcomplex > &Hamiltonian_1,
               Eigen::SparseMatrix< Eigen::dcomplex > &Hamiltonian_2,
               const int &div_x, const double &x_min, const double &x_max, const double &dt, const double &mass,
               const double &width, const std::complex<double> &strength, const QVector< std::complex<double> > &otherParameter,
               std::complex<double> (*potential)(const double &, const double &, const std::complex<double> &, const QVector< std::complex<double> > &))
{
    double dx = (x_max - x_min)/div_x;
    const std::complex<double> constant1 = ( std::complex<double>(0.0, 1.0) * dt ) / ( 2.0 * hc );
    const std::complex<double> constant2 = -constant1 * hc * hc / ( 2.0 * mass * dx * dx);
    Eigen::SparseMatrix<Eigen::dcomplex> I(div_x, div_x); /* unit matrix */
    Eigen::VectorXcd V(div_x);
    QVector<Eigen::Triplet<Eigen::dcomplex>> triplets; /* elements list of sparse matrix */

    /* calc potential vector */
    for(int ix = 0; ix < div_x; ++ix){
        V(ix) = potential(x_min + dx * ix, width, strength, otherParameter);
    }

    /* calc hamiltonian matrix */

/*    triplets.emplaceBack(0, 0, 0);
    triplets.emplaceBack(0, 1, 0);
    triplets.emplaceBack(div_x - 1, div_x - 2, 0);
    triplets.emplaceBack(div_x - 1, div_x - 1, 0);
    for(int row = 1; row < div_x - 1; ++row){
        for(int clm = row - 1; clm <= row + 1 ; ++clm){
            // hamiltonian matrix is tridiagonal matrix
            if(( clm < 0 ) || ( clm >= div_x )){
                continue;
            }
            triplets.emplace_back(row, clm,
                                  constant2 * (double)(delta(row, clm - 1) - 2.0 * delta(row, clm) + delta(row, clm + 1) )
                                      + constant1 * V(row) * (double)delta(row, clm) );
        }
    }
*/
    for(int row = 0; row < div_x; ++row){
        for(int clm = row - 1; clm <= row + 1 ; ++clm){
            // hamiltonian matrix is tridiagonal matrix
            if(( clm < 0 ) || ( clm >= div_x )){
                continue;
            }
            triplets.emplace_back(row, clm,
                                  constant2 * (double)(delta(row, clm - 1) - 2.0 * delta(row, clm) + delta(row, clm + 1) )
                                + constant1 * V(row) * (double)delta(row, clm) );
        }
    }
    /* calc coefficient matrix */
    Hamiltonian_1.resize(div_x, div_x);
    Hamiltonian_2.resize(div_x, div_x);
    I.setIdentity();
    Hamiltonian_2.setFromTriplets(triplets.begin(), triplets.end());
    Hamiltonian_1 = I - Hamiltonian_2;
    Hamiltonian_2 = I + Hamiltonian_2;

    /* set solver */
    solver_prev.compute(Hamiltonian_1);
    solver_next.compute(Hamiltonian_2);
    if(( solver_prev.info() != Eigen::Success ) || ( solver_next.info() != Eigen::Success )){
        return false;
    }
    return true;
}


/**
 * @brief time evolution function
 *
 * @param[in,out] wave_function in:current wave functino, out:wave function after time evolution.
 * @param[in] solver
 * @param[in] Hamiltonian coefficient matrix for current wave function.
 *
 * @return true:success, false:failure
 */
bool time_evolution(Eigen::VectorXcd &wave_function,
                      const Eigen::SparseQR< Eigen::SparseMatrix<Eigen::dcomplex>, Eigen::COLAMDOrdering<int> > &solver,
                      const Eigen::SparseMatrix< Eigen::dcomplex > &Hamiltonian)
{
    Eigen::VectorXcd temp = Hamiltonian * wave_function;
    wave_function = solver.solve(temp);
    if(solver.info() != Eigen::Success){
        return false;
    }
    return true;
}

/**
 * @brief set Gaussian wave packet
 *
 * @param[out] wave_function
 * @param[in] x0 center of wave packet
 * @param[in] p momentum of wave packet
 * @param[in] width width of wave packet
 * @param[in] div_x stride of x
 * @param[in] x_min min of xAxis
 * @param[in] x_max max of xAxis
 */
void set_Gaussian_wave_packet(Eigen::VectorXcd &wave_function, const double &x0, const std::complex<double> &p, const double &width,
                                    const int &div_x, const double &x_min, const double &x_max)
{
    double x;
    double dx = (x_max - x_min) / div_x;
    double normalized_constant = pow(2.0 * pi * width * width, -0.25);

    wave_function.resize(div_x);
    for(int ix = 0; ix < div_x; ++ix){
        x = x_min + dx * ix;
        wave_function(ix) = normalized_constant * exp( - pow((x - x0) / (2.0 * width), 2) + std::complex<double>(0.0, 1.0) * p * (x - x0) / hc );
    }
}

/**
 * @brief calc relative mass
 *
 * @param[in] my_n
 * @param[in] target_n
 * @param[in] target_p
 *
 * @return relative mass
 */
double ret_mass(const int &my_n, const int &target_n, const int &target_p)
{
    return mn * my_n * ( mn * target_n + mp * target_p ) / ( mn * ( my_n + target_n) + mp * target_p );
}

/**
 * @brief calc momentum
 * @param[in] E
 * @param[in] mass
 *
 * @return momentum
 */
std::complex<double> ret_p(const std::complex<double> &E, const double &mass)
{
    return std::sqrt(2.0 * mass * E);
}
