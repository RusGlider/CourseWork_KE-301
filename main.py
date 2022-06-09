from PyQt5.QtWidgets import*
from PyQt5.uic import loadUi
from matplotlib.backends.backend_qt5agg import (NavigationToolbar2QT as NavigationToolbar)
import numpy as np
import widget_plot
import ode_system

class MainWidget(QMainWindow):
    solver = ode_system.ODESolver()
    def __init__(self):
        QMainWindow.__init__(self)

        loadUi("windowdesign.ui", self)

        self.setWindowTitle("Модель Центрального Регуляторного Контура Системы Морфогенеза Механорецепторов Дрозофилы")
        self.Button_run.clicked.connect(self.run)
        self.Button_default.clicked.connect(self.set_defaults)
        #self.pushButton_generate_random_signal.clicked.connect(self.update_graph)
        self.addToolBar(NavigationToolbar(self.widget_plot.canvas, self))
        self.set_interface_values()

    def run(self):
        self.set_system_values()
        time = float(self.lineEdit_t.text())
        t,Y = self.solver.make_experiment(time)
        self.draw_graph(t,Y)

    def check_inputs(self):
        pass


    def set_system_values(self):

        new_concs, new_coeffs, new_mutations, new_sigmoids, new_s = ([] for i in range(5))

        new_concs.append(float(self.lineEdit_D.text()))
        new_concs.append(float(self.lineEdit_E.text()))
        new_concs.append(float(self.lineEdit_G.text()))
        new_concs.append(float(self.lineEdit_S.text()))
        new_concs.append(float(self.lineEdit_U.text()))
        new_concs.append(float(self.lineEdit_C.text()))
        new_concs.append(float(self.lineEdit_L.text()))

        new_coeffs.append(float(self.lineEdit_kx.text()))
        new_coeffs.append(float(self.lineEdit_ky.text()))
        new_coeffs.append(float(self.lineEdit_kz.text()))
        new_coeffs.append(float(self.lineEdit_ku.text()))
        new_coeffs.append(float(self.lineEdit_kw.text()))
        new_coeffs.append(float(self.lineEdit_kp.text()))

        new_mutations.append(float(self.lineEdit_mx.text()))
        new_mutations.append(float(self.lineEdit_my.text()))
        new_mutations.append(float(self.lineEdit_mz.text()))
        new_mutations.append(float(self.lineEdit_mu.text()))
        new_mutations.append(float(self.lineEdit_mw.text()))
        new_mutations.append(float(self.lineEdit_mp.text()))
        new_mutations.append(float(self.lineEdit_m7.text()))

        new_sigmoids.append(float(self.lineEdit_n1.text()))
        new_sigmoids.append(float(self.lineEdit_d1.text()))
        new_sigmoids.append(float(self.lineEdit_n3.text()))
        new_sigmoids.append(float(self.lineEdit_d3.text()))
        new_sigmoids.append(float(self.lineEdit_n5.text()))
        new_sigmoids.append(float(self.lineEdit_d5.text()))

        new_s.append(float(self.lineEdit_a3.text()))
        new_s.append(float(self.lineEdit_b3.text()))
        new_s.append(float(self.lineEdit_c3.text()))
        new_s.append(float(self.lineEdit_a4.text()))
        new_s.append(float(self.lineEdit_b4.text()))
        new_s.append(float(self.lineEdit_c4.text()))
        new_s.append(float(self.lineEdit_a5.text()))
        new_s.append(float(self.lineEdit_b5.text()))
        new_s.append(float(self.lineEdit_c5.text()))
        new_s.append(float(self.lineEdit_a6.text()))
        new_s.append(float(self.lineEdit_b6.text()))
        new_s.append(float(self.lineEdit_c6.text()))

        new_tau = float(self.lineEdit_tau.text())
        new_time = float(self.lineEdit_t.text())

        self.solver.update_parameters(new_concs, new_coeffs, new_mutations, new_sigmoids, new_s, new_tau, new_time)


    def set_interface_values(self):
        self.lineEdit_D.setText(str(self.solver.D))
        self.lineEdit_G.setText(str(self.solver.G))
        self.lineEdit_E.setText(str(self.solver.E))
        self.lineEdit_S.setText(str(self.solver.S))
        self.lineEdit_U.setText(str(self.solver.U))
        self.lineEdit_L.setText(str(self.solver.L))
        self.lineEdit_C.setText(str(self.solver.C))

        self.lineEdit_kx.setText(str(self.solver.kx))
        self.lineEdit_ky.setText(str(self.solver.ky))
        self.lineEdit_kz.setText(str(self.solver.kz))
        self.lineEdit_ku.setText(str(self.solver.ku))
        self.lineEdit_kw.setText(str(self.solver.kw))
        self.lineEdit_kp.setText(str(self.solver.kp))

        self.lineEdit_mx.setText(str(self.solver.mx))
        self.lineEdit_my.setText(str(self.solver.my))
        self.lineEdit_mz.setText(str(self.solver.mz))
        self.lineEdit_mu.setText(str(self.solver.mu))
        self.lineEdit_mw.setText(str(self.solver.mw))
        self.lineEdit_mp.setText(str(self.solver.mp))
        self.lineEdit_m7.setText(str(self.solver.m7))

        self.lineEdit_n1.setText(str(self.solver.n1))
        self.lineEdit_d1.setText(str(self.solver.d1))
        self.lineEdit_n3.setText(str(self.solver.n3))
        self.lineEdit_d3.setText(str(self.solver.d3))
        self.lineEdit_n5.setText(str(self.solver.n5))
        self.lineEdit_d5.setText(str(self.solver.d5))

        self.lineEdit_a3.setText(str(self.solver.a3))
        self.lineEdit_b3.setText(str(self.solver.b3))
        self.lineEdit_c3.setText(str(self.solver.c3))
        self.lineEdit_a4.setText(str(self.solver.a4))
        self.lineEdit_b4.setText(str(self.solver.b4))
        self.lineEdit_c4.setText(str(self.solver.c4))
        self.lineEdit_a5.setText(str(self.solver.a5))
        self.lineEdit_b5.setText(str(self.solver.b5))
        self.lineEdit_c5.setText(str(self.solver.c5))
        self.lineEdit_a6.setText(str(self.solver.a6))
        self.lineEdit_b6.setText(str(self.solver.b6))
        self.lineEdit_c6.setText(str(self.solver.c6))

        self.lineEdit_t.setText(str(self.solver.default_t))
        self.lineEdit_tau.setText(str(self.solver.tau))

    def draw_graph(self, t, Y):
        self.widget_plot.canvas.axes.clear()
        handles = []

        if (self.checkBox_ASC.isChecked()):
            self.widget_plot.canvas.axes.plot(t, Y[:, 0], 'r-', linewidth=2.0, label="ASC(t)")
            handles.append('ASC(t)')
        if (self.checkBox_Hairy.isChecked()):
            self.widget_plot.canvas.axes.plot(t, Y[:, 1], 'g--', linewidth=2.0, label="Hairy(t)")
            handles.append('Hairy(t)')
        if (self.checkBox_SENS.isChecked()):
            self.widget_plot.canvas.axes.plot(t, Y[:, 2], 'b-.', linewidth=2.0, label="SENS(t)")
            handles.append('SENS(t)')
        if (self.checkBox_SCRT.isChecked()):
            self.widget_plot.canvas.axes.plot(t, Y[:, 3], 'm-', linewidth=2.0, label="SCRT(t)")
            handles.append('SCRT(t)')
        if (self.checkBox_CHN.isChecked()):
            self.widget_plot.canvas.axes.plot(t, Y[:, 4], 'y--', linewidth=2.0, label="CHN(t)")
            handles.append('CHN(t)')
        if (self.checkBox_PHYL.isChecked()):
            self.widget_plot.canvas.axes.plot(t, Y[:, 5], 'c-.', linewidth=2.0, label="PHYL(t)")
            handles.append('PHYL(t)')

        self.widget_plot.canvas.axes.set_xlabel("time")
        self.widget_plot.canvas.axes.set_ylabel("concentration")

        self.widget_plot.canvas.axes.legend(handles, loc='upper right')
        self.widget_plot.canvas.axes.set_title('Концентрации белков')
        self.widget_plot.canvas.draw()

    def set_defaults(self):
        self.solver.set_default_parameters()
        self.set_interface_values()

def main():
    app = QApplication([])
    window = MainWidget()
    window.show()
    app.exec_()

if __name__ == "__main__":
    main()