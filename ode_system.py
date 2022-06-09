import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

class ODESolver:
    """
    # Концентрации белков
    D = 2.05   #   DA
    E = 1   #   EMC
    G = 1   #   GRO
    S = 5.6   #   SINA
    U = 1.17   #   Ub
    C = 1   #   SCRT ?
    L = 1   #

    # Скорости деградации соответствующих белков
    mx = 1.77
    my = 1
    mz = 1
    mu = 1
    mw = 1
    mp = 1
    m7 = 1

    # Коэффициенты мутации
    kx = 1
    ky = 1
    kz = 1
    ku = 1
    kw = 1
    kp = 1

    # Коэффициенты сигмоидных функций
    n1 = 2
    d1 = 7.46

    n3 = 2
    d3 = 2.77

    n5 = 2
    d5 = 1.24

    # Коэффициенты функции S
    a3 = 3.61
    b3 = 4.96
    c3 = 1.35

    a4 = 4.43
    b4 = 6.09
    c4 = 1.66

    a5 = 8.09
    b5 = 11.13
    c5 = 3.03

    a6 = 2.67
    b6 = 3.67
    c6 = 1

    tau = 12
    default_t = 17
    """
    def __init__(self):
        self.set_default_parameters()

    def set_default_parameters(self):
        # Концентрации белков
        self.D = 2.05  # DA
        self.E = 1  # EMC
        self.G = 1  # GRO
        self.S = 5.6  # SINA
        self.U = 1.17  # Ub
        self.C = 1  # SCRT ?
        self.L = 1  #

        # Скорости деградации соответствующих белков
        self.mx = 1.77
        self.my = 1
        self.mz = 1
        self.mu = 1
        self.mw = 1
        self.mp = 1
        self.m7 = 1

        # Коэффициенты мутации
        self.kx = 1
        self.ky = 1
        self.kz = 1
        self.ku = 1
        self.kw = 1
        self.kp = 1

        # Коэффициенты сигмоидных функций
        self.n1 = 2
        self.d1 = 7.46

        self.n3 = 2
        self.d3 = 2.77

        self.n5 = 2
        self.d5 = 1.24

        # Коэффициенты функции S
        self.a3 = 3.61
        self.b3 = 4.96
        self.c3 = 1.35

        self.a4 = 4.43
        self.b4 = 6.09
        self.c4 = 1.66

        self.a5 = 8.09
        self.b5 = 11.13
        self.c5 = 3.03

        self.a6 = 2.67
        self.b6 = 3.67
        self.c6 = 1

        self.tau = 12
        self.default_t = 17

    def update_parameters(self, concs, coeffs, mutations, sigmoids, s, tau, time):
        self.update_protein_concentrations(concs)
        self.update_mutation_coeffs(coeffs)
        self.update_protein_degradation_speeds(mutations)
        self.update_sigmoid_params(sigmoids)
        self.update_s_params(s)
        self.update_tau(tau)
        self.default_t = time

    def update_protein_concentrations(self, concs):
        self.D = concs[0]   #new_D
        self.E = concs[1]   #new_E
        self.G = concs[2]   #new_G
        self.S = concs[3]   #new_S
        self.U = concs[4]   #new_U
        self.C = concs[5]   #new_C
        self.L = concs[6]   #new_L

    def update_protein_degradation_speeds(self, mutations): # new_mx, new_my, new_mz, new_mu, new_mw, new_mp, new_m7):
        self.mx = mutations[0]   #new_mx
        self.my = mutations[1]   #new_my
        self.mz = mutations[2]   #new_mz
        self.mu = mutations[3]   #new_mu
        self.mw = mutations[4]   #new_mw
        self.mp = mutations[5]   #new_mp
        self.m7 = mutations[6]   #new_m7

    def update_mutation_coeffs(self, coeffs): #new_kx, new_ky, new_kz, new_ku, new_kw, new_kp):
        self.kx = coeffs[0]   #new_kx
        self.ky = coeffs[1]   #new_ky
        self.kz = coeffs[2]   #new_kz
        self.ku = coeffs[3]   #new_ku
        self.kw = coeffs[4]   #new_kw
        self.kp = coeffs[5]   #new_kp

    def update_sigmoid_params(self, sigmoids): #new_n1, new_d1, new_n3, new_d3, new_n5, new_d5):
        self.n1 = sigmoids[0]   #new_n1
        self.d1 = sigmoids[1]   #new_d1
        self.n3 = sigmoids[2]   #new_n3
        self.d3 = sigmoids[3]   #new_d3
        self.n5 = sigmoids[4]   #new_n5
        self.d5 = sigmoids[5]   #new_d5

    def update_s_params(self, s): # new_a3, new_b3, new_c3,new_a4 ,new_b4 ,new_c4, new_a5, new_b5, new_c5, new_a6, new_b6, new_c6):
        self.a3 = s[0]  #new_a3
        self.b3 = s[1]  #new_b3
        self.c3 = s[2]  #new_c3
        self.a4 = s[3]  #new_a4
        self.b4 = s[4]  #new_b4
        self.c4 = s[5]  #new_c4
        self.a5 = s[6]  #new_a5
        self.b5 = s[7]  #new_b5
        self.c5 = s[8]  #new_c5
        self.a6 = s[9]  #new_a6
        self.b6 = s[10]  #new_b6
        self.c6 = s[11]  #new_c6

    def update_tau(self, new_tau):
        #global tau
        self.tau = new_tau
    """
    def System(self, X, t):
        dxdt = ( kx*(sigma1(D*X[0])+sigma3(X[2])+sigma5(X[4])) / ((1+G*X[1])*(1+E*X[0])) ) - (1+X[5]*(t-tau)*U*S)*mx*X[0]
        dydt = (ky*(C/(d1+(X[3]**m7))))-my*X[1]
        dzdt = kz*s3(D*X[0])-mz*X[2]
        dudt = ku*s4(D*X[0])-mu*X[3]
        dwdt = kw*s5(D*X[0])-mw*X[4]
        dpdt = kp*((s6(D*X[0])*heaviside(t-tau)*((t-tau)**2))/(L + heaviside(t-tau)*((t-tau)**2))) - mp*X[5]
        return [dxdt,dydt,dzdt,dudt,dwdt,dpdt]
    """

    def System(self, X, t):
        dxdt = (self.kx * (self.sigma1(self.D * X[0]) + self.sigma3(X[2]) + self.sigma5(X[4])) /
                ((1 + self.G * X[1]) * (1 + self.E * X[0]))) - (1 + X[5] * (t - self.tau) * self.U * self.S) * self.mx * X[0]
        dydt = (self.ky * (self.C / (self.d1 + (X[3] ** self.m7)))) - self.my * X[1]
        dzdt = self.kz * self.s3(self.D * X[0]) - self.mz * X[2]
        dudt = self.ku * self.s4(self.D * X[0]) - self.mu * X[3]
        dwdt = self.kw * self.s5(self.D * X[0]) - self.mw * X[4]
        dpdt = self.kp * ((self.s6(self.D * X[0]) * self.heaviside(t - self.tau) * ((t - self.tau) ** 2)) / (
                    self.L + self.heaviside(t - self.tau) * ((t - self.tau) ** 2))) - self.mp * X[5]
        return [dxdt, dydt, dzdt, dudt, dwdt, dpdt]

    def ExperimentXYZ(self):
        t = np.linspace(0, 17, 1000)
        X0 = [0.56, 1.59, 0.15, 0, 0, 0]
        Y = odeint(self.System, X0, t)

        plt.plot(t, Y[:, 0], 'r-', linewidth=2.0, label="x(t)")
        plt.plot(t, Y[:, 1], 'g--', linewidth=2.0, label="y(t)")
        plt.plot(t, Y[:, 2], 'b-.', linewidth=2.0, label="z(t)")
        plt.xlabel("t")
        plt.ylabel("x(t), y(t), z(t)")
        plt.legend()
        plt.grid()
        plt.show()

    def ExperimentUWP(self):
        t = np.linspace(0, 17, 1000)
        X0 = [0.56, 1.59, 0.15, 0, 0, 0]
        Y = odeint(self.System, X0, t)

        plt.plot(t, Y[:, 3], 'r-', linewidth=2.0, label="u(t)")
        plt.plot(t, Y[:, 4], 'g--', linewidth=2.0, label="w(t)")
        plt.plot(t, Y[:, 5], 'b-.', linewidth=2.0, label="p(t)")
        plt.xlabel("t")
        plt.ylabel("u(t), w(t), p(t)")
        plt.legend()
        plt.grid()
        plt.show()

    def make_experiment(self, time = 15):
        t = np.linspace(0, time, 1000)
        X0 = [0.56, 1.59, 0.15, 0, 0, 0]
        Y = odeint(self.System, X0, t)
        return [t,Y]


    def sigma1(self, q):
        return self.d1 * (q ** self.n1) / (1 + (q ** self.n1))

    def sigma3(self, q):
        return self.d3 * (q ** self.n3) / (1 + (q ** self.n3))

    def sigma5(self, q):
        return self.d5 * (q ** self.n5) / (1 + (q ** self.n5))

    #def sigma(q, n, d):
        #return d*(q**n) / (1 + (q**n))

    def s3(self, q):
        return (self.a3*np.exp((q-self.b3)/self.c3)) / (1+np.exp((q-self.b3)/self.c3))

    def s4(self, q):
        return (self.a4*np.exp((q-self.b4)/self.c4)) / (1+np.exp((q-self.b4)/self.c4))

    def s5(self, q):
        return (self.a5*np.exp((q-self.b5)/self.c5)) / (1+np.exp((q-self.b5)/self.c5))

    def s6(self, q):
        return (self.a6*np.exp((q-self.b6)/self.c6)) / (1+np.exp((q-self.b6)/self.c6))
    """
    def s(q,a,b,c):
        return (a*math.exp((q-b)/c)) / (1+math.exp((q-b)/c))
    """
    def heaviside(self, t):
        return 0 if t < 0 else 1
