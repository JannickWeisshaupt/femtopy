import numpy as np
import matplotlib.pyplot as plt


class GaussBeam:
    def __init__(self, lambda0, v1, v2=0, mode='w&R', Ep=1, fwhm_t=1):
        self.lambda_0 = lambda0
        self.Ep = Ep
        self.fwhm_t = fwhm_t
        self.k = 2 * np.pi / lambda0
        self.P = Ep / (1.06 * fwhm_t)

        if mode.lower() == 'w&r':
            win = v1
            Rin = v2
            self.q0 = (1 / Rin - 1j * 2 / (self.k * win ** 2)) ** (-1)
        elif mode.lower() == "q":
            self.q0 = v1
        elif mode.lower() == 'w&z':
            win = v1
            zin = v2
            self.q0 = zin + 1j * (self.k * win ** 2 / 4 + np.sqrt((self.k * win ** 2 / 4) ** 2 - zin ** 2))

    def propagate(self, mges):
        self.q0 = (mges[0, 0] * self.q0 + mges[0, 1]) / (mges[1, 0] * self.q0 + mges[1, 1])

    def w(self):
        return np.sqrt(2 / (self.k * np.abs((self.q0 ** -1).imag)))

    def R(self):
        return 1 / (self.q0 ** -1).real

    def z(self):
        return self.q0.real

    def z0(self):
        return self.q0.imag

    def w0(self):
        return self.w() / np.sqrt(1 + (np.pi * self.w() ** 2 / (self.lambda_0 * self.R())) ** 2)

    def phi(self):
        return 2 * self.lambda_0 / (np.pi * self.w0())

    def Imax(self):
        return 2 * self.P / (np.pi * self.w() ** 2)

    def Imax0(self):
        return 2 * self.P / (np.pi * self.w0() ** 2)

    def __repr__(self):
        return 'Actual Position: z = {0:2.3f} m.  w = {1:2.2e} m.    R = {2:1.2e} m.  Imax = {3:1.2e} W/cm^2\nFocus:           z0= {4:2.2f} m.  w0= {5:1.2e} m.  phi = {6:1.2e} rad.   Imax0 = {7:1.2e} W/cm^2  \n'.format(
            *[self.z(), self.w(), self.R(), self.Imax() * 1e-4, self.z0(), self.w0(), self.phi(), self.Imax0() * 1e-4])


def lens(f):
    return np.array([[1, 0], [-1 / f, 1]])


def dist(d):
    return np.array([[1, d], [0, 1]])


if __name__ == '__main__':
    b1 = GaussBeam(400e-9, 0.05, 100000000)
    b2 = GaussBeam(266e-9, 0.05, 100000000)

