import abc
import numpy as np
from scipy.interpolate import interp1d


class Profile(abc.ABC):
    """Represent a density profile, allowing it to be sampled from"""

    def __init__(self, r_cutoff, **kwargs):
        self.r_cutoff = r_cutoff

    @abc.abstractmethod
    def density(self, r):
        """The density profile as a function of r"""
        pass

    def cdf(self, ndim):
        """
        Returns the cumulative distribution function of the density - taking into account how this varies
        in n dimensions
        """
        r = np.linspace(0, self.r_cutoff, 1000)
        rho = self.density(r)
        # dn = rho * dr * r**(N-1)
        cdf = np.cumsum(rho * r ** (ndim - 1))
        if cdf[0] != 0:
            cdf -= cdf[0]
        cdf /= cdf.max()
        return r, cdf

    def sample(self, n, ndim):
        r, cdf = self.cdf(ndim)

        interpolator = interp1d(cdf, r)

        # Sampling from an arbitrary distribution is equivalent to sampling from a uniform distribution then
        # reverse mapping the output number in the range 0-1 to the CDF of the desired distribution
        samples = interpolator(np.random.random(n))

        return samples


class ConstantProfile(Profile):
    """Represent a constant density"""

    def density(self, r):
        return np.ones_like(r)


class ExponentialProfile(Profile):
    """Represent an exponential distribution"""

    def __init__(self, r_cutoff, h):
        self.r_cutoff = r_cutoff
        self.h = h

    def density(self, r):
        return np.exp(-r / self.h)


class PlummerProfile(Profile):
    """Represent an Plummer distribution"""

    def __init__(self, r_cutoff, a):
        self.r_cutoff = r_cutoff
        self.a = a

    def density(self, r):
        return 1 / self.a**3 * (1 + r**2 / self.a**2) ** (-5 / 2)


class NFWProfile(Profile):
    """Represent a Navarro Frenk White profile"""

    def __init__(self, r_cutoff, a):
        self.r_cutoff = r_cutoff
        self.a = a

    def density(self, r):
        return 1 / (r / self.a) / (1 + r / self.a) ** 2


class HernquistProfile(Profile):
    """Represent a Hernquist profile"""

    def __init__(self, r_cutoff, a):
        self.r_cutoff = r_cutoff
        self.a = a

    def density(self, r):
        return self.a / r / (r + self.a) ** 3


class IsothermalProfile(Profile):
    """Represernt an isothermal profile"""

    def __init__(self, r_cutoff, a):
        self.r_cutoff = r_cutoff
        self.a = a

    def density(self, r):
        return 1 / (1 + (r / self.a) ** 2)
