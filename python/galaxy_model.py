import abc
import numpy as np

from state import State
from profiles import Profile


class _ThetaProfile(Profile):
    """For sampling the polar angle correctly in 3D"""

    def __init__(self):
        self.r_cutoff = np.pi

    def density(self, r):
        return np.sin(r)


class GalaxyComponent(abc.ABC):
    """Represent the component of a galaxy (e.g., bulge, disk etc)"""

    ndim = None

    def __init__(self, n_objs, profile: Profile, mass=None):
        self.n_objs = n_objs
        self.mass = self.n_objs if mass is None else mass
        self.state = State.empty(n_objs, ndim=self.ndim)
        self.profile = profile

    @abc.abstractmethod
    def star_positions(self):
        """Method to compute the positions of particles"""
        pass


class Isotropic2D(GalaxyComponent):
    """Galaxy profile that only depends on r (2-dimensional)"""

    ndim = 2

    def star_positions(self):
        r = self.profile.sample(self.n_objs, self.ndim)
        phi = np.random.uniform(0, 2 * np.pi, self.n_objs)

        x = r * np.cos(phi)
        y = r * np.sin(phi)

        return np.asarray([x, y]).T


class Isotropic3D(GalaxyComponent):
    """Galaxy profile that only depends on r (3-dimensional)"""

    ndim = 3

    def star_positions(self):

        r = self.profile.sample(self.n_objs, self.ndim)
        # theta must be sampled such that P(theta) = sin(theta)
        theta = _ThetaProfile().sample(self.n_objs, ndim=1)
        phi = np.random.uniform(0, 2 * np.pi, self.n_objs)

        x = r * np.sin(theta) * np.cos(phi)
        y = r * np.sin(theta) * np.sin(phi)
        z = r * np.cos(theta)

        return np.asarray([x, y, z]).T


class Disk2D(GalaxyComponent):
    """2D only disk component"""

    ndim = 2

    def __init__(self, n_objs, disk_profile, mass=None):
        self.n_objs = n_objs
        self.disk_profile = disk_profile
        self.mass = self.n_objs if mass is None else mass

    def star_positions(self):
        r = self.disk_profile.sample(self.n_objs, self.ndim)
        phi = np.random.uniform(0, 2 * np.pi, self.n_objs)

        x = r * np.cos(phi)
        y = r * np.sin(phi)

        return np.asarray([x, y]).T


class Disk3D(GalaxyComponent):
    """3D disk component with separate profiles for the surface density and height density"""

    ndim = 3

    def __init__(self, n_objs, disk_profile, height_profile, mass=None):
        self.n_objs = n_objs
        self.disk_profile = disk_profile
        self.height_profile = height_profile
        self.mass = self.n_objs if mass is None else mass

    def star_positions(self):
        r = self.disk_profile.sample(self.n_objs, ndim=2)
        z = self.height_profile.sample(self.n_objs, ndim=1) * np.random.choice((-1, 1), self.n_objs)
        phi = np.random.uniform(0, 2 * np.pi, self.n_objs)

        x = r * np.cos(phi)
        y = r * np.sin(phi)

        return np.asarray([x, y, z]).T
