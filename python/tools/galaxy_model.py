import abc
import numpy as np
import matplotlib.pyplot as plt

from .state import State
from .profiles import _Profile


class _ThetaProfile(_Profile):
    """For sampling the polar angle correctly in 3D"""

    def __init__(self):
        self.r_cutoff = np.pi

    def density(self, r):
        return np.sin(r)


class _GalaxyComponent(abc.ABC):
    """Represent the component of a galaxy (e.g., bulge, disk etc)"""

    ndim = None

    def __init__(self, n_objs, profile: _Profile, mass=None, velocity_dispersion=0.1, tag=0):
        self.n_objs = n_objs
        self.mass = self.n_objs if mass is None else mass
        self.state = State.empty(n_objs, ndim=self.ndim)
        self.profile = profile
        self.velocity_dispersion = velocity_dispersion
        self.tag = tag

    @abc.abstractmethod
    def star_positions(self):
        """Method to compute the positions of particles"""
        pass

    @abc.abstractmethod
    def star_velocities(self, positions):
        """Compute velocity directions (e.g. unit velocities according to positions)"""


class Isotropic2D(_GalaxyComponent):
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

        flip = np.random.choice([-1, 1], self.n_objs)

        u *= flip
        v *= flip

        dispersion = np.random.normal(0, 1, (self.n_objs, self.ndim))
        vel = np.asarray([u, v]).T
        vel += self.velocity_dispersion * dispersion

        return vel


class Isotropic3D(_GalaxyComponent):
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
        xhat[:, 0] = 1
        v1 = np.cross(positions, xhat)
        v2 = np.cross(positions, v1)

        # Normalise
        v1 = v1 / np.linalg.norm(v1, axis=-1, keepdims=True)
        v2 = v2 / np.linalg.norm(v2, axis=-1, keepdims=True)

        # Randomise the orientations
        theta = np.random.uniform(0, 2 * np.pi, self.n_objs)
        vel = v1 * np.cos(theta[:, np.newaxis]) + v2 * np.sin(theta[:, np.newaxis])

        dispersion = np.random.normal(0, 1, (self.n_objs, self.ndim))
        vel += self.velocity_dispersion * dispersion
        return vel


class Disk2D(_GalaxyComponent):
    """2D only disk component"""

    ndim = 2

    def __init__(self, n_objs, disk_profile, mass=None, velocity_dispersion=0.1, tag=0):
        self.n_objs = n_objs
        self.disk_profile = disk_profile
        self.mass = self.n_objs if mass is None else mass
        self.velocity_dispersion = velocity_dispersion
        self.tag = tag

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

        dispersion = np.random.normal(0, 1, (self.n_objs, self.ndim))
        vel = np.asarray([u, v]).T
        vel += self.velocity_dispersion * dispersion

        return vel


class Disk3D(_GalaxyComponent):
    """3D disk component with separate profiles for the surface density and height density"""

    ndim = 3

    def __init__(self, n_objs, disk_profile, height_profile, mass=None, velocity_dispersion=0.1, tag=0):
        self.n_objs = n_objs
        self.disk_profile = disk_profile
        self.height_profile = height_profile
        self.mass = self.n_objs if mass is None else mass
        self.velocity_dispersion = velocity_dispersion
        self.tag = tag

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

        dispersion = np.random.normal(0, 1, (self.n_objs, self.ndim))
        vel = np.asarray([u, v, w]).T
        vel += self.velocity_dispersion * dispersion

        return vel


class Galaxy:
    def __init__(self, state):
        self.state = state

    @classmethod
    def from_components(cls, components_list):
        """Combine input components and compute velocities"""
        positions_list = []
        velocities_list = []
        masses_list = []
        tags_list = []

        # concatenate components
        for component in components_list:
            pos = component.star_positions()
            vel = component.star_velocities(pos)
            masses = np.ones(component.n_objs) * np.float64(component.mass) / np.float64(component.n_objs)
            tags = np.full(component.n_objs, component.tag, dtype=np.int32)

            positions_list.append(pos)
            velocities_list.append(vel)
            masses_list.append(masses)
            tags_list.append(tags)

        pos = np.vstack(positions_list)
        vel = np.vstack(velocities_list)
        mass = np.hstack(masses_list)
        tag = np.hstack(tags_list)

        # create state
        state = State.empty(len(pos))
        state.pos[:, :] = pos[:, :]
        state.vel[:, :] = vel[:, :]
        state.mass[:] = mass[:]
        state.tag[:] = tag[:]

        # Compute velocity magnitudes and apply
        r = np.sqrt((state.pos**2).sum(-1))

        sorted_inds = np.argsort(r)
        m_cumulative = np.cumsum(state.mass[sorted_inds])
        # Circular velocity = sqrt(GM(r)/r)
        vc = np.sqrt(m_cumulative / r[sorted_inds])

        vels = np.zeros_like(state.mass)
        vels[sorted_inds] = vc

        state.vel *= vels[:, np.newaxis]

        # remove mean (mass weighted) velocity from system
        v_mean = np.sum(state.vel*state.mass[:, np.newaxis], axis=0) / state.mass.sum()
        state.vel -= v_mean[np.newaxis, :]

        com = np.sum(state.pos*state.mass[:, np.newaxis], axis=0)/state.mass.sum()
        state.pos -= com[np.newaxis, :]

        return Galaxy(state)

    def normalise_masses_for_time(self):
        """Adjust masses so t=1 represents a crossing time"""
        # t_c = sqrt(r^3/M) -> M = r**3/1**2
        m = self.state.mass.sum()
        r = np.std(self.state.pos, axis=0).sum()

        m_req = r**3

        self.state.mass *= m_req / m

        # also need to modify velocities
        # v^2 = GM/R -> velocity needs to be scaled by sqrt(m)
        self.state.vel *= np.sqrt(m_req/m)

        print(f"Old mass = {m}, new_mass = {self.state.mass.sum()}")

    def write(self, filename=None, normalise_time=True):
        if normalise_time:
            self.normalise_masses_for_time()
        print("Writing galaxy to file")
        self.state.write(filename)

    def shift(self, pos, vel):
        """Shift whole galaxy"""

        self.state.pos += pos[np.newaxis, :]
        self.state.vel += vel[np.newaxis, :]

    def rotate(self, xz_rot, xy_rot):
        """Apply a rotation first about the xz plane, then about the xy plane"""

        xz = np.radians(xz_rot)
        xy = np.radians(xy_rot)

        r1 = np.asarray(
            [
                [np.cos(xz), 0, -np.sin(xz)],
                [0, 1, 0],
                [np.sin(xz), 0, np.cos(xz)],
            ]
        )

        r2 = np.asarray(
            [
                [np.cos(xy), -np.sin(xy), 0],
                [np.sin(xy), np.cos(xy), 0],
                [0, 0, 1],
            ]
        )

        r = r1 @ r2

        self.state.pos = self.state.pos @ r
        self.state.vel = self.state.vel @ r

    def set_experiment_name(self, name):
        self.state.name = name

    def __add__(self, other):
        return Galaxy(self.state + other.state)
