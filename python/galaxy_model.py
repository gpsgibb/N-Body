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

    @abc.abstractmethod
    def star_velocities(self, positions):
        """Compute velocity directions (e.g. unit velocities according to positions)"""


class Isotropic2D(GalaxyComponent):
    """Galaxy profile that only depends on r (2-dimensional)"""

    ndim = 2

    def star_positions(self):
        r = self.profile.sample(self.n_objs, self.ndim)
        phi = np.random.uniform(0, 2 * np.pi, self.n_objs)

        x = r * np.cos(phi)
        y = r * np.sin(phi)

        return np.asarray([x, y]).T

    def star_velocities(self, positions):
        # We want velocities orthogonal to r
        x, y = positions.T
        theta = np.arctan2(y, x)
        u = -np.sin(theta)
        v = np.cos(theta)

        flip = np.random.choice([-1,1], self.n_objs)

        u *= flip
        v *= flip

        return np.asarray([u, v]).T


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


    def star_velocities(self, positions):
        # Cross r with x unit vector to get a vector orthogonal to r and x
        # (so a vector on the tangent plane of the sphere)
        xhat = np.zeros((self.n_objs, 3))
        xhat[:,0] = 1
        v1 = np.cross(positions, xhat)
        v2 = np.cross(positions, v1)

        # Normalise
        v1 = v1/np.linalg.norm(v1, axis=-1, keepdims=True)
        v2 = v2/np.linalg.norm(v2, axis=-1, keepdims=True)

        # Randomise the orientations
        theta = np.random.uniform(0, 2*np.pi, self.n_objs)
        vel = v1*np.cos(theta[:, np.newaxis]) + v2*np.sin(theta[:, np.newaxis])

        return vel


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

    def star_velocities(self, positions):
        x, y = positions.T
        theta = np.arctan2(y, x)

        u = -np.sin(theta)
        v = np.cos(theta)

        return np.asarray([u, v]).T


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

    def star_velocities(self, positions):
        x, y, z = positions.T
        theta = np.arctan2(y, x)

        u = -np.sin(theta)
        v = np.cos(theta)
        w = np.zeros_like(u)

        return np.asarray([u, v, w]).T
