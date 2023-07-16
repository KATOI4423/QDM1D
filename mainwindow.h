#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include "SparseCore"
#include "SparseQR"

QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

private slots:
    void on_bWell_clicked();
    void on_bWS_clicked();
    void on_lineReE_editingFinished();
    void on_lineImE_editingFinished();
    void on_lineR0_editingFinished();
    void on_lineReV_editingFinished();
    void on_lineImV_editingFinished();
    void on_lineA_editingFinished();
    void on_lineProton_editingFinished();
    void on_lineNeutron_editingFinished();
    void on_lineX0_editingFinished();
    void on_lineWaveWidth_editingFinished();
    void on_barTime_valueChanged(int value);
    void on_lineTmax_editingFinished();
    void on_lineDt_editingFinished();
    void on_lineXmin_editingFinished();
    void on_lineXmax_editingFinished();
    void on_lineYmin_editingFinished();
    void on_lineYmax_editingFinished();
    void on_b2RG_clicked();

    void on_lineReV1_editingFinished();

    void on_lineImV1_editingFinished();

    void on_lineMu0_editingFinished();

    void on_lineMu1_editingFinished();

private:
    Ui::MainWindow *ui;

    /* Axis */
    int div_x = 1001;
    double x_min;
    double x_max;
    double y_min;
    double y_max;
    QVector<double> x;

    /* slider */
    int t_min = 0;
    int t_max;
    int t_now = 0;
    double dt;

    /* for Quantum Dynamics Mechanics */
    Eigen::SparseQR< Eigen::SparseMatrix<Eigen::dcomplex>, Eigen::COLAMDOrdering<int> > solver_prev; /* solver for dt < 0 */
    Eigen::SparseQR< Eigen::SparseMatrix<Eigen::dcomplex>, Eigen::COLAMDOrdering<int> > solver_next; /* solver for dt > 0 */
    Eigen::SparseMatrix< Eigen::dcomplex > Hamiltonian_1; /* H for psi(t) */
    Eigen::SparseMatrix< Eigen::dcomplex > Hamiltonian_2; /* H for psi(t + dt) */
    Eigen::VectorXcd psi; /* wave function */
    double x0;
    double wave_width;
    double expX;
    double expPsi;

    /* nuclear */
    int target_n;
    int target_p;
    int my_n;
    double mu;
    std::complex<double> E;
    std::complex<double> p;

    /* potential */
    double r0;
    double potential_width;
    std::complex<double> potensial_strength;
    QVector<std::complex<double>> potential_otherParameter;
    std::complex<double> (*potential)(const double &x, const double &width, const std::complex<double> &strength, const QVector< std::complex<double> > &otherParameter); /* function pointer of potential */

    /* method */
    void addGraph(void);
    void initialize(void);
};
#endif // MAINWINDOW_H
