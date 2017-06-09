import numpy as np
import matplotlib.pyplot as plt


class Sellmeyer:
    """Class Sellmeyer is a Sellmeyer type refractive index representation
    and can be used to calculate optical properties from Sellmeyer or cauchy coefficients"""

    def __init__(self, n0_sq=1.0, resonances_sell=(), resonances_cauchy=(), linear=0.0, quadratic=0.0):
        self.n0 = n0_sq
        self.resonances = resonances_sell
        self.resonances_cauchy = resonances_cauchy
        self.linear = linear
        self.quadratic = quadratic

        self.c = 299792458

    def refractive_index(self,lambda0):
        return np.sqrt(self.refractive_index_squared(lambda0))

    def refractive_index_squared(self, lambda0):
        nout_sq = self.n0
        for el in self.resonances:
            nout_sq += el[0] * lambda0 ** 2 / (lambda0 ** 2 - el[1])
        for el_c in self.resonances_cauchy:
            nout_sq += el_c[0]/ (lambda0 ** 2 - el_c[1])
        nout_sq += self.quadratic * lambda0 ** 2
        nout_sq += self.linear * lambda0
        return nout_sq

    def del_n_sq2_del_lambda2(self, lambda0):
        res = 0
        for el in self.resonances:
            res += 2 * el[0]*el[1]*(el[1] +3*lambda0**2 ) / (lambda0 ** 2 - el[1]) ** 3
        for el_c in self.resonances_cauchy:
            res += 2 * el_c[0] * (el_c[1]+3*lambda0**2) / (lambda0 ** 2 - el_c[1]) ** 3
        res += self.quadratic*2
        return res


    def del_n2_del_lambda2(self, lambda0):
        return 0.5*self.del_n_sq2_del_lambda2(lambda0)/self.refractive_index(lambda0)-0.25*self.refractive_index_squared(lambda0)**(-3/2)*(self.del_n_sq_del_lambda(lambda0))**2

    def del_n_sq_del_lambda(self, lambda0):
        res = 0
        for el in self.resonances:
            res += -2 * el[0] * el[1] * lambda0 / (lambda0 ** 2 - el[1]) ** 2
        for el_c in self.resonances_cauchy:
            res += -2 * el_c[0] * lambda0 / (lambda0 ** 2 - el_c[1]) ** 2
        res += 2 * self.quadratic * lambda0
        res += self.linear
        return res

    def del_n_del_lambda(self, lambda0):
        return 0.5 * self.del_n_sq_del_lambda(lambda0) / self.refractive_index(lambda0)

    def group_velocity_index(self,lambda0):
        return self.refractive_index(lambda0)-lambda0*self.del_n_del_lambda(lambda0)

    def group_velocity(self,lambda0):
        return self.c/self.group_velocity_index(lambda0)

    def D(self, lambda0):
        return -lambda0*1e6/self.c*self.del_n2_del_lambda2(lambda0)

    def gvd(self, lambda0,unit='SI'):
        if unit not in ['SI','fs^2/mm']:
            raise Exception('unit must be SI or fs^2/mm')
        res = -(lambda0*1e-6)**2/(2*np.pi*self.c)*self.D(lambda0)
        if unit == 'fs^2/mm':
            res = res*1e27
        return res

if __name__ == "__main__":

    lambda0 = np.linspace(0.15, 1, 1000)

    kbbf_o = Sellmeyer(1,[[1.169725,0.0062400]],quadratic=-0.009904)
    kbbf_e = Sellmeyer(1, [[0.956611, 0.0061926]], quadratic=-0.027849)

    n_kbbf_o = kbbf_o.refractive_index(lambda0)
    n_kbbf_e = kbbf_e.refractive_index(lambda0)

    plt.figure(1)
    plt.plot(lambda0,n_kbbf_e)
    plt.plot(lambda0,n_kbbf_o)

    plt.figure(2)
    plt.plot(lambda0,kbbf_e.gvd(lambda0))
    plt.plot(lambda0,kbbf_o.gvd(lambda0))

    plt.figure(3)
    plt.plot(lambda0,kbbf_e.group_velocity_index(lambda0))
    plt.plot(lambda0,kbbf_o.group_velocity_index(lambda0))

    # #### BBO
    # bbo_o = Sellmeyer(n0_sq=2.7405, resonances_cauchy=[[0.0184, 0.0179]], quadratic=-0.0155)
    #
    #
    # plt.figure(1)
    # # plt.plot(lambda0,n_kbbf_o,lambda0,n_kbbf_e)
    # plt.plot(lambda0, bbo_o.refractive_index(lambda0))
    #
    # plt.ylim(1.625, 1.85)
    #
    # plt.figure(2)
    # plt.plot(lambda0,bbo_o.group_velocity_index(lambda0))
    #
    # print(bbo_o.group_velocity_index(lambda0.min()))
    #
    #
    # plt.figure(3)
    # plt.plot(lambda0,bbo_o.gvd(lambda0))


    # print(bbo_o.D(lambda0.min()))