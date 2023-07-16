#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "quantum_dynamics_1d.h"
#include <QMessageBox>

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    addGraph();
    this->potential = potential_well;
    initialize();
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::addGraph(void)
{
    double str;
    /* グラフ追加 */
    ui->widget->addGraph(); /* graph(0) = Re(psi) */
    ui->widget->addGraph(); /* graph(1) = Im(psi) */
    ui->widget->addGraph(ui->widget->xAxis, ui->widget->yAxis2); /* graph(2) = Re(V) */
    ui->widget->addGraph(ui->widget->xAxis, ui->widget->yAxis2); /* graph(3) = Im(V) */
    ui->widget->addGraph(); /* graph(4) = (psi, psi) */

    ui->widget->yAxis2->setVisible(true);
    ui->widget->graph(0)->setName("Re(ψ)");
        ui->widget->graph(1)->setName("Im(ψ)");
        ui->widget->graph(2)->setName("Re(V)");
    ui->widget->graph(3)->setName("Im(V)");
    ui->widget->graph(4)->setName("|ψ|^2");

        /* グラフの色 */
        QPen pen;
    pen.setWidth(2);
    pen.setColor(QColor(255,0,0));
    ui->widget->graph(0)->setPen(pen);
    pen.setColor(QColor(0,0,255));
    ui->widget->graph(1)->setPen(pen);
    pen.setColor(QColor(0,0,0));
    ui->widget->graph(2)->setPen(pen);
    pen.setColor(QColor(100,100,100));
    pen.setStyle(Qt::DashLine);
    ui->widget->graph(3)->setPen(pen);
    pen.setColor(QColor(100, 100, 255));
    pen.setStyle(Qt::DotLine);
    ui->widget->graph(4)->setPen(pen);

    /* ラベル名 */
    ui->widget->xAxis->setLabel("x [fm]");
    ui->widget->yAxis->setLabel("Wave function");
    ui->widget->yAxis2->setLabel("Potential [MeV]");

    /* 凡例 */
    ui->widget->legend->setVisible(true);
    ui->widget->legend->setBrush(QBrush(Qt::white));
    ui->widget->legend->setFillOrder(QCPLegend::foColumnsFirst);
    ui->widget->axisRect()->insetLayout()->setInsetAlignment(0, Qt::AlignTop | Qt::AlignLeft);
    ui->widget->plotLayout()->setMargins(QMargins(0, 0, 0, 0));
    for(int i=0; i<5; ++i){
        QCPAbstractLegendItem* item = new QCPPlottableLegendItem(ui->widget->legend, ui->widget->graph(i));
        ui->widget->legend->addItem(item);
        delete item;
    }
}

void MainWindow::initialize(void)
{
    double dx;
    double temp = 0;
    QVector<std::complex<double>> V;
    QVector<double> RePsi, ImPsi, Dens, ReV, ImV;

    /* set parameter */
    this->t_now = 0;
    ui->barTime->setValue(0);
    ui->gTime->setTitle("time = 0");
    this->dt = ui->lineDt->text().toDouble();

    this->x_max = ui->lineXmax->text().toDouble();
    this->x_min = ui->lineXmin->text().toDouble();
    this->y_max = ui->lineYmax->text().toDouble();
    this->y_min = ui->lineYmin->text().toDouble();
    dx = (this->x_max - this->x_min) / this->div_x;

    this->target_n = ui->lineNeutron->text().toInt();
    this->target_p = ui->lineProton->text().toInt();
    this->mu = ret_mass(this->my_n, this->target_n, this->target_p);
    this->E = std::complex<double>(ui->lineReE->text().toDouble(), ui->lineImE->text().toDouble());
    this->p = ret_p(E, this->mu);
    this->x0 = ui->lineX0->text().toDouble();
    this->wave_width = ui->lineWaveWidth->text().toDouble();

    this->r0 = ui->lineR0->text().toDouble();
    if(this->potential_otherParameter.size() == 0){
        this->potential_otherParameter.resize(4);
    }
    this->potential_otherParameter[0] = ui->lineA->text().toDouble();
    this->potential_otherParameter[1] = std::complex<double>(ui->lineReV1->text().toDouble(), ui->lineImV1->text().toDouble());
    this->potential_otherParameter[2] = ui->lineMu0->text().toDouble();
    this->potential_otherParameter[3] = ui->lineMu1->text().toDouble();
    this->potential_width = this->r0 * cbrt(ui->lineNeutron->text().toInt() + ui->lineProton->text().toInt());
    this->potensial_strength = std::complex<double>(ui->lineReV->text().toDouble(), ui->lineImV->text().toDouble());

    /* set vector size */
    x.resize(this->div_x);
    V.resize(this->div_x);
    RePsi.resize(this->div_x);
    ImPsi.resize(this->div_x);
    Dens.resize(this->div_x);
    ReV.resize(this->div_x);
    ImV.resize(this->div_x);

    /* set init wave packet */
    set_Gaussian_wave_packet(this->psi, this->x0, this->p, this->wave_width, this->div_x, this->x_min, this->x_max);

    /* calc */
    for(int ix = 0; ix < this->div_x; ++ix){
        this->x[ix] = this->x_min + dx * ix;
        V[ix] = this->potential(x[ix], this->potential_width, this->potensial_strength, this->potential_otherParameter);
        RePsi[ix] = this->psi(ix).real();
        ImPsi[ix] = this->psi(ix).imag();
        Dens[ix] = pow(abs(this->psi(ix)), 2);
        ReV[ix] = V[ix].real();
        ImV[ix] = V[ix].imag();
        if(abs(V[ix].real()) > temp){
            temp = abs(V[ix].real());
        }
        if(abs(V[ix].imag()) > temp){
            temp = abs(V[ix].imag());
        }
    }

    /* set axis */
    ui->widget->xAxis->setRange(this->x_min, this->x_max);
    ui->widget->yAxis->setRange(this->y_min, this->y_max);
    ui->widget->yAxis2->setRange(-temp, temp);

    /* set data */
    ui->widget->graph(0)->setData(this->x, RePsi);
    ui->widget->graph(1)->setData(this->x, ImPsi);
    ui->widget->graph(2)->setData(this->x, ReV);
    ui->widget->graph(3)->setData(this->x, ImV);
    ui->widget->graph(4)->setData(this->x, Dens);
    ui->widget->replot();

    /* set Infomation */
    this->expX = this->x0;
    this->expPsi = 1.0;
    ui->lexpX->setText("<x> [fm] = " + QString::number(this->expX, 'g', 16));
    ui->lexpPsi->setText("<ψ|ψ> = 1.0");

    ui->lhc->setText("hc = " + QString::number(hc, 'g', 16) + " [Mev・s / c^2]");
    ui->lmp->setText("m_p = " + QString::number(mp, 'g', 16) + " [MeV / c^2]");
    ui->lmn->setText("m_n = " + QString::number(mn, 'g', 16) + " [MeV / c^2]");
    ui->lDx->setText("dx = " + QString::number(dx, 'g', 16) + " [fm]");

    /* set solver */
    set_Solver(this->solver_prev, this->solver_next, this->Hamiltonian_1, this->Hamiltonian_2, this->div_x, this->x_min, this->x_max, this->dt, this->mu,
               this->potential_width, this->potensial_strength, this->potential_otherParameter, this->potential);
}

/* on_function */

void MainWindow::on_bWell_clicked()
{
    ui->bWS->setChecked(false);
    ui->b2RG->setChecked(false);
    ui->lA->setEnabled(false);
    ui->lineA->setEnabled(false);
    ui->lReV1->setEnabled(false);
    ui->lineReV1->setEnabled(false);
    ui->lImV1->setEnabled(false);
    ui->lineImV1->setEnabled(false);
    ui->lMu0->setEnabled(false);
    ui->lineMu0->setEnabled(false);
    ui->lMu1->setEnabled(false);
    ui->lineMu1->setEnabled(false);

    this->potential = potential_well;
    initialize();
}

void MainWindow::on_bWS_clicked()
{
    ui->bWell->setChecked(false);
    ui->b2RG->setChecked(false);
    ui->lA->setEnabled(true);
    ui->lineA->setEnabled(true);
    ui->lReV1->setEnabled(false);
    ui->lineReV1->setEnabled(false);
    ui->lImV1->setEnabled(false);
    ui->lineImV1->setEnabled(false);
    ui->lMu0->setEnabled(false);
    ui->lineMu0->setEnabled(false);
    ui->lMu1->setEnabled(false);
    ui->lineMu1->setEnabled(false);

    this->potential = woods_saxon;
    initialize();
}


void MainWindow::on_b2RG_clicked()
{
    ui->bWell->setChecked(false);
    ui->bWS->setChecked(false);
    ui->lA->setEnabled(false);
    ui->lineA->setEnabled(false);
    ui->lReV1->setEnabled(true);
    ui->lineReV1->setEnabled(true);
    ui->lImV1->setEnabled(true);
    ui->lineImV1->setEnabled(true);
    ui->lMu0->setEnabled(true);
    ui->lineMu0->setEnabled(true);
    ui->lMu1->setEnabled(true);
    ui->lineMu1->setEnabled(true);

    this->potential = two_range_gausian;
    ui->lineReV->setText(QString::number(this->potential_otherParameter[1].real() * -2.0, 'g', 6));
    ui->lineImV->setText(QString::number(this->potential_otherParameter[1].imag() * -2.0, 'g', 6));
    initialize();
}

void MainWindow::on_lineReE_editingFinished()
{
    bool check;
    double value = ui->lineReE->text().toDouble(&check);
    if(check == false){
        ui->lineReE->setText(QString::number(this->E.real(), 'g', 16));
        return;
    }
    if(value == this->E.real()){
        return;
    }
    initialize();
}

void MainWindow::on_lineImE_editingFinished()
{
    bool check;
    double value = ui->lineImE->text().toDouble(&check);
    if(check == false){
        ui->lineImE->setText(QString::number(this->E.imag(), 'g', 16));
        return;
    }
    if(value == this->E.imag()){
        return;
    }
    initialize();
}

void MainWindow::on_lineR0_editingFinished()
{
    bool check;
    double value = ui->lineR0->text().toDouble(&check);
    if(( check == false ) || ( value <= 0 )){
        ui->lineR0->setText(QString::number(this->r0, 'g', 16));
        return;
    }
    if(value == this->r0){
        return;
    }
    initialize();
}

void MainWindow::on_lineReV_editingFinished()
{
    bool check;
    double value = ui->lineReV->text().toDouble(&check);
    if(check == false){
        ui->lineReV->setText(QString::number(this->potensial_strength.real(), 'g', 16));
        return;
    }
    if(value == this->potensial_strength.real()){
        return;
    }
    initialize();
}

void MainWindow::on_lineImV_editingFinished()
{
    bool check;
    double value = ui->lineImV->text().toDouble(&check);
    if(check == false){
        ui->lineImV->setText(QString::number(this->potensial_strength.imag(), 'g', 16));
        return;
    }
    if(value == this->potensial_strength.imag()){
        return;
    }
    initialize();
}

void MainWindow::on_lineA_editingFinished()
{
    bool check;
    double value = ui->lineA->text().toDouble(&check);
    if(( check == false ) || ( value <= 0 )){
        ui->lineA->setText(QString::number(this->potential_otherParameter[0].real(), 'g', 16));
        return;
    }
    if(value == this->potential_otherParameter[0].real()){
        return;
    }
    initialize();
}

void MainWindow::on_lineReV1_editingFinished()
{
    bool check;
    double value = ui->lineReV1->text().toDouble(&check);
    if(check == false){
        ui->lineReV1->setText(QString::number(this->potential_otherParameter[1].real(), 'g', 16));
        return;
    }
    if(value == this->potential_otherParameter[1].real()){
        return;
    }
    initialize();
}

void MainWindow::on_lineImV1_editingFinished()
{
    bool check;
    double value = ui->lineImV1->text().toDouble(&check);
    if(check == false){
        ui->lineImV1->setText(QString::number(this->potential_otherParameter[1].imag(), 'g', 16));
        return;
    }
    if(value == this->potential_otherParameter[1].imag()){
        return;
    }
    initialize();
}

void MainWindow::on_lineMu0_editingFinished()
{
    bool check;
    double value = ui->lineMu0->text().toDouble(&check);
    if(check == false){
        ui->lineMu0->setText(QString::number(this->potential_otherParameter[2].real(), 'g', 16));
        return;
    }
    if(value == this->potential_otherParameter[2].real()){
        return;
    }
    initialize();
}

void MainWindow::on_lineMu1_editingFinished()
{
    bool check;
    double value = ui->lineMu1->text().toDouble(&check);
    if(check == false){
        ui->lineMu1->setText(QString::number(this->potential_otherParameter[3].real(), 'g', 16));
        return;
    }
    if(value == this->potential_otherParameter[3].real()){
        return;
    }
    initialize();
}


void MainWindow::on_lineProton_editingFinished()
{
    bool check;
    int value = ui->lineProton->text().toInt(&check);
    if(( check == false ) || ( value <= 0 ) ||
        ( ( ui->lineProton->text().toInt() == 0 ) && ( ui->lineNeutron->text().toInt() == 0 ) )){
        ui->lineProton->setText(QString::number(this->target_p, 10));
        return;
    }
    if(value == this->target_p){
        return;
    }
    initialize();
}

void MainWindow::on_lineNeutron_editingFinished()
{
    bool check;
    int value = ui->lineNeutron->text().toInt(&check);
    if(( check == false ) || ( value <= 0 ) ||
        ( ( ui->lineProton->text().toInt() == 0 ) && ( ui->lineNeutron->text().toInt() == 0 ) )){
        ui->lineNeutron->setText(QString::number(this->target_n, 10));
        return;
    }
    if(value == this->target_n){
        return;
    }
    initialize();
}

void MainWindow::on_lineX0_editingFinished()
{
    bool check;
    double value = ui->lineX0->text().toDouble(&check);
    if(( check == false ) || ( value <= 0 ) || (value > this->x_max )){
        ui->lineX0->setText(QString::number(this->x0, 'g', 16));
        return;
    }
    if(value == this->x0){
        return;
    }
    initialize();
}

void MainWindow::on_lineWaveWidth_editingFinished()
{
    bool check;
    double value = ui->lineWaveWidth->text().toDouble(&check);
    if(( check == false ) || ( value <= 0 )){
        ui->lineWaveWidth->setText(QString::number(this->wave_width, 'g', 16));
        return;
    }
    if(value == this->potential_width){
        return;
    }
    initialize();
}

void MainWindow::on_lineTmax_editingFinished()
{
    bool check;
    int value = ui->lineTmax->text().toInt(&check);
    if(( check == false ) || ( value <= this->t_min )){
        ui->lineTmax->setText(QString::number(this->t_max, 10));
        return;
    }
    this->t_max = ui->lineTmax->text().toInt();
    ui->barTime->setRange(this->t_min, this->t_max);
}

void MainWindow::on_lineDt_editingFinished()
{
    bool check;
    double value = ui->lineDt->text().toDouble(&check);
    if(( check == false ) || ( value <= 0 )){
        ui->lineDt->setText(QString::number(this->dt, 'g', 16));
        return;
    }
    if(value == this->dt){
        return;
    }
    initialize();
}

void MainWindow::on_lineXmin_editingFinished()
{
    bool check;
    double value = ui->lineXmin->text().toDouble(&check);
    if(( check == false ) || ( value >= this->x_max )){
        ui->lineXmin->setText(QString::number(this->x_min, 'g', 6));
        return;
    }
    if(value == this->x_min){
        return;
    }
    if(this->x0 <= value){
        this->x0 = this->x_min + (this->x_max - this->x_min) * 0.3;
        ui->lineX0->setText(QString::number(this->x0, 'g', 6));
    }
    initialize();
}

void MainWindow::on_lineXmax_editingFinished()
{
    bool check;
    double value = ui->lineXmax->text().toDouble(&check);
    if(( check == false ) || ( value <= this->x_min )){
        ui->lineXmax->setText(QString::number(this->x_max, 'g', 16));
        return;
    }
    if(value == this->x_max){
        return;
    }
    if(this->x0 >= value){
        this->x0 = this->x_min + (this->x_max - this->x_min) * 0.3;
        ui->lineX0->setText(QString::number(this->x0, 'g', 6));
    }
    initialize();
}

void MainWindow::on_lineYmin_editingFinished()
{
    bool check;
    double value = ui->lineYmin->text().toDouble(&check);
    if(( check == false ) || ( value >= this->y_max )){
        ui->lineYmin->setText(QString::number(this->y_min, 'g', 16));
        return;
    }
    if(value == this->y_min){
        return;
    }
    initialize();
}

void MainWindow::on_lineYmax_editingFinished()
{
    bool check;
    double value = ui->lineYmax->text().toDouble(&check);
    if(( check == false ) || ( value <= this->y_min )){
        ui->lineYmax->setText(QString::number(this->y_max, 'g', 16));
        return;
    }
    if(value == this->y_max){
        return;
    }
    initialize();
}

void MainWindow::on_barTime_valueChanged(int value)
{
    int operation;
    double dx = ( this->x_max - this->x_min ) / this->div_x;
    Eigen::SparseQR< Eigen::SparseMatrix<Eigen::dcomplex>, Eigen::COLAMDOrdering<int> > *solver;
    Eigen::SparseMatrix< Eigen::dcomplex > *H;
    QVector<double> RePsi(this->div_x), ImPsi(this->div_x), Dens(this->div_x);

    if(value > this->t_now){
        solver = &(this->solver_next);
        H = &(this->Hamiltonian_1);
        operation = +1;
    }else if(value < this->t_now){
        solver = &(this->solver_prev);
        H = &(this->Hamiltonian_2);
        operation = -1;
    }else{
        return;
    }

    do{
        this->expX = 0.0;
        this->expPsi = 0.0;

        if(time_evolution(this->psi, *solver, *H) == false){
            QMessageBox::critical(this, "ERROR OF SOLVER", "Solver cannot calculation this time evolution.");
            initialize();
            return;
        }
        for(int ix = 0; ix < this->div_x; ++ix){
            RePsi[ix] = this->psi(ix).real();
            ImPsi[ix] = this->psi(ix).imag();
            Dens[ix] = pow(abs(this->psi(ix)), 2);
            if(ix == 0) continue;

            this->expX += ( ( this->x_min + dx * ix ) * Dens[ix] + ( this->x_min + dx * ( ix - 1) ) * Dens[ix - 1] ) * 0.5 * dx;
            this->expPsi += ( Dens[ix] + Dens[ix - 1] ) * 0.5  * dx;
        }
        ui->widget->graph(0)->setData(this->x, RePsi);
        ui->widget->graph(1)->setData(this->x, ImPsi);
        ui->widget->graph(4)->setData(this->x, Dens);
        ui->widget->replot();
        this->t_now += operation;
        ui->lexpX->setText("<x> [fm] = " + QString::number(this->expX, 'g', 16));
        ui->lexpPsi->setText("<ψ|ψ> = " + QString::number(this->expPsi, 'g', 16));
        ui->gTime->setTitle("time = " + QString::number(this->t_now, 10));
    }while(this->t_now != value);
}

