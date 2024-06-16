import numpy as np
from scipy import interpolate
from scipy.integrate import quad

MINVAL = 1e-7


def func_zeros(a, b, func, delta=1e-7):
    if func(a) * func(b) > 0:
        print(a, func(a), b, func(b))
        raise ValueError("ends have the same sign")
    while abs(b - a) > delta:
        y1 = func(a)
        y2 = func(0.5 * (a + b))
        y3 = func(b)

        if y1 * y2 < 0:
            b = 0.5 * (a + b)
        else:
            a = 0.5 * (a + b)
    return (a + b) * 0.5


class FlatBrush:
    def __init__(
        self, N: int, sigma: float, chi: float, eta: float = 1, num_bins: int = 100
    ) -> None:
        self.N = N
        self.sigma = sigma
        self.chi = chi
        self.eta = eta
        self.H = float(N) * 0.5
        self.zt = 0.0
        self.Hmax = float(N) - MINVAL
        self.num_bins = num_bins
        _ = func_zeros(MINVAL, self.Hmax, self.target_H)
        self.phi = np.hstack([self.phi, np.array([0.0, 0.0])])
        self.z = np.hstack([self.z, np.array([self.H + 0.5*self.dz, self.H+1.5*self.dz])])

    def target(self, phi):
        return (
            3 / 2 * (np.pi * self.eta / (2 * self.N)) ** 2 * (self.H**2 - self.zt**2)
            + np.log(1.0 - phi)
            + 2 * self.chi * phi
        )

    def phi_on_z(self, z):
        self.zt = z
        return func_zeros(0.0, 1.0 - MINVAL, self.target)

    def target_H(self, H):
        self.H = H
        self.z = np.linspace(0.0, H, self.num_bins)[1:]
        self.dz = H / (self.num_bins - 1)
        self.phi = np.array([self.phi_on_z(i) for i in self.z - self.dz / 2])
        return np.sum(self.phi) * self.dz - self.N * self.sigma

    @property
    def hight(self):
        return self.H
    
    
class TwoBrush:    
    def __init__(self, N: int, sigma: float, chi: float, eta: float = 1, num_bins: int = 100) -> None:
        self.N = N
        self.sigma = sigma
        self.chi = chi
        self.eta = eta
        self.num_bins = num_bins
        fb = FlatBrush(N, sigma, chi, eta=eta, num_bins=num_bins)
        self.phi = interpolate.interp1d(fb.z, fb.phi, fill_value = "extrapolate")
        self.d = np.linspace(N*sigma+fb.dz, fb.H, num_bins)
        self.H = fb.H
        
    def phi_midplace(self, d):
        return self.phi(d) + quad(self.phi, d, self.N, args=())[0] / d
    
    def overlap_zone(self, a, d):
        return a * self.N**(2/3) * (2*d)**(-1/3) 
    
    def pressure(self, d):
        phi = self.phi_midplace(d)
        return  -phi - np.log(1 - phi) - self.chi * phi**2
    
    def overlap_integral(self, a, d):
        return self.phi_midplace(d)**2 * self.overlap_zone(a, d)
        
        
            
