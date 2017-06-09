import numpy as np
import matplotlib.pyplot as plt


class Sellmeyer:
    """Class Sellmeyer is a Sellmeyer type refractive index representation
    and can be used to calculate optical properties from Sellmeyer or cauchy coefficients

    The refractive index is given by the formula:

    n^2 = n0_sq + Sum_i A_i lambda^2/(lambda^2-B_i) + Sum_j C_j/(lambda^2-D_j) +  quadratic lambda^2 + linear * lambda

    where A_i and B_i are the sellmeyer coefficients
    C_j and D_j are the Cauchy coefficients
    n0_sq, linear and quadratic are polynomial coeffiecients

    The input arguments are:
    n0_sq=1.0               constant term
    linear=0.0              linear term
    quadratic=0.0           quadratic term
    resonances_sell=()      List of two elements list for sellmeyer coefficients (Example: [[1,0.5],[2,0.1]]
    resonances_cauchy=()    List of two elements list for cauchy coefficients (Example: [[1,0.5],[2,0.1]]

    Example usage:

    bbo_o = Sellmeyer(n0_sq=2.7405, resonances_cauchy=[[0.0184, 0.0179]], quadratic=-0.0155)
    lambda0 = np.linspace(0.2,1,1000) # Lambda is typically in micrometers
    n_bbo_o = bbo_o.refractive_index(lambda0)
    """

    def __init__(self, n0_sq=1.0, resonances_sell=(), resonances_cauchy=(), linear=0.0, quadratic=0.0):
        self.n0 = n0_sq
        self.resonances = resonances_sell
        self.resonances_cauchy = resonances_cauchy
        self.linear = linear
        self.quadratic = quadratic

        self.c = 299792458

    def refractive_index(self, lambda0):
        """"Function to calculate the refractive index """

        lambda0 = self._convert_to_array(lambda0)
        return np.sqrt(self.refractive_index_squared(lambda0))

    def refractive_index_squared(self, lambda0):
        lambda0 = self._convert_to_array(lambda0)
        nout_sq = self.n0
        for el in self.resonances:
            nout_sq += el[0] * lambda0 ** 2 / (lambda0 ** 2 - el[1])
        for el_c in self.resonances_cauchy:
            nout_sq += el_c[0] / (lambda0 ** 2 - el_c[1])
        nout_sq += self.quadratic * lambda0 ** 2
        nout_sq += self.linear * lambda0
        return nout_sq

    def del_n_sq2_del_lambda2(self, lambda0):
        lambda0 = self._convert_to_array(lambda0)
        res = 0
        for el in self.resonances:
            res += 2 * el[0] * el[1] * (el[1] + 3 * lambda0 ** 2) / (lambda0 ** 2 - el[1]) ** 3
        for el_c in self.resonances_cauchy:
            res += 2 * el_c[0] * (el_c[1] + 3 * lambda0 ** 2) / (lambda0 ** 2 - el_c[1]) ** 3
        res += self.quadratic * 2
        return res

    def del_n2_del_lambda2(self, lambda0):
        lambda0 = self._convert_to_array(lambda0)
        return 0.5 * self.del_n_sq2_del_lambda2(lambda0) / self.refractive_index(
            lambda0) - 0.25 * self.refractive_index_squared(lambda0) ** (-3 / 2) * (self.del_n_sq_del_lambda(
            lambda0)) ** 2

    def del_n_sq_del_lambda(self, lambda0):
        lambda0 = self._convert_to_array(lambda0)
        res = 0
        for el in self.resonances:
            res += -2 * el[0] * el[1] * lambda0 / (lambda0 ** 2 - el[1]) ** 2
        for el_c in self.resonances_cauchy:
            res += -2 * el_c[0] * lambda0 / (lambda0 ** 2 - el_c[1]) ** 2
        res += 2 * self.quadratic * lambda0
        res += self.linear
        return res

    def del_n_del_lambda(self, lambda0):
        lambda0 = self._convert_to_array(lambda0)
        return 0.5 * self.del_n_sq_del_lambda(lambda0) / self.refractive_index(lambda0)

    def group_velocity_index(self, lambda0):
        lambda0 = self._convert_to_array(lambda0)
        return self.refractive_index(lambda0) - lambda0 * self.del_n_del_lambda(lambda0)

    def group_velocity(self, lambda0):
        lambda0 = self._convert_to_array(lambda0)
        return self.c / self.group_velocity_index(lambda0)

    def D(self, lambda0):
        lambda0 = self._convert_to_array(lambda0)
        return -lambda0 * 1e6 / self.c * self.del_n2_del_lambda2(lambda0)

    def gvd(self, lambda0, unit='SI'):
        lambda0 = self._convert_to_array(lambda0)
        if unit not in ['SI', 'fs^2/mm']:
            raise Exception('unit must be SI or fs^2/mm')
        res = -(lambda0 * 1e-6) ** 2 / (2 * np.pi * self.c) * self.D(lambda0)
        if unit == 'fs^2/mm':
            res = res * 1e27
        return res

    @staticmethod
    def _convert_to_array(input_var):
        if type(input_var) in [tuple,list]:
            output = np.array(input_var, dtype=np.float)
        elif type(input_var) in [np.ndarray, float, int]:
            output = input_var
        else:
            raise ValueError('Input wavelength must be float, list of floats, or np.array')
        return output

if __name__ == "__main__":
    lambda_test = np.linspace(0.15, 1, 1000)

    kbbf_o = Sellmeyer(1, [[1.169725, 0.0062400]], quadratic=-0.009904)
    kbbf_e = Sellmeyer(1, [[0.956611, 0.0061926]], quadratic=-0.027849)

    n_kbbf_o = kbbf_o.refractive_index(lambda_test)
    n_kbbf_e = kbbf_e.refractive_index(lambda_test)

    plt.figure(1)
    plt.plot(lambda_test, n_kbbf_e)
    plt.plot(lambda_test, n_kbbf_o)

    plt.figure(2)
    plt.plot(lambda_test, kbbf_e.gvd(lambda_test))
    plt.plot(lambda_test, kbbf_o.gvd(lambda_test))

    plt.figure(3)
    plt.plot(lambda_test, kbbf_e.group_velocity_index(lambda_test))
    plt.plot(lambda_test, kbbf_o.group_velocity_index(lambda_test))

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
