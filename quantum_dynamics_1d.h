/**
 * @code quantum_dynamics_1d.h
 * @brief calc 1 Dimension Quantum Dynamics with Crank-Nicolson method
 */

#ifndef QUANTUM_DYNAMICS_1D_H
#define QUANTUM_DYNAMICS_1D_H

#include <QVector>
#include "SparseCore"
#include "SparseQR"

constexpr double pi = 4.0 * atan(1.0);
constexpr double hc = 197.3271 /* planck constant [MeV・s / c^2] */;
constexpr double mp = 938.27208816; /* mass for proton [MeV / c^2]*/
constexpr double mn = 939.5654133; /* mass for neutron [MeV / c^2]*/

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
                                    const std::complex<double> &strength, const QVector<std::complex<double>> &otherParameter);


/**
 * @brief woods_saxon potential
 *
 * @param[in] x
 * @param[in] width nuclear radius.
 * @param[in] strength strength of potential. strength > 0 means it is hight. strength < 0 means it is depth.
 * @param[in] otherParameter [0]:a
 *
 * @return optential at x.
 */
std::complex<double> woods_saxon(const double &x, const double &width,
                                 const std::complex<double> &strength, const QVector<std::complex<double>> &otherParameter);

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
                                       const std::complex<double> &strength, const QVector<std::complex<double>> &otherParameter);

/**
 * @brief set Solver for 1 Direction Quantum Dynamics
 *
 * @param[out] solver_prev solver for t - dt
 * @param[out] solver_next solver for t + dt
 * @param[out] Hamiltonian_1 coefficient matrix for psi(t)
 * @param[out] Hamiltonian_2 coefficient matrix for psi(t + dt)
 * @param[in] div_x stride of x
 * @param[in] x_min min of xAxis
 * @param[in] x_max max of xAxis
 * @param[in] dt stride of time
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
               std::complex<double> (*potential)(const double &, const double &, const std::complex<double> &, const QVector< std::complex<double> > &));

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
                    const Eigen::SparseMatrix< Eigen::dcomplex > &Hamiltonian);

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
                              const int &div_x, const double &x_min, const double &x_max);

/**
 * @brief calc relative mass
 *
 * @param[in] my_n
 * @param[in] target_n
 * @param[in] target_p
 *
 * @return relative mass
 */
double ret_mass(const int &my_n, const int &target_n, const int &target_p);

/**
 * @brief calc momentum
 * @param[in] E
 * @param[in] mass
 *
 * @return momentum
 */
std::complex<double> ret_p(const std::complex<double> &E, const double &mass);

#endif // QUANTUM_DYNAMICS_1D_H
