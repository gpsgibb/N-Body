from dataclasses import dataclass

import numpy as np
import h5py


@dataclass
class State:
    mass: np.ndarray
    pos: np.ndarray
    vel: np.ndarray
    acc: np.ndarray
    pe: np.ndarray
    ke: np.ndarray
    density: np.ndarray
    pressure: np.ndarray
    smoothing_length: np.ndarray
    tag: np.ndarray
    ndim: int
    iteration: int = 0
    time: float = 0.0
    name: str = "default"

    @classmethod
    def read(cls, filename):
        with h5py.File(filename, "r") as h:

            time = h.attrs["time"]
            iteration = h.attrs["iteration"]
            ndim = h.attrs["NDIM"]
            name = h.attrs["name"]

            pos = h["pos"][:]
            vel = h["vel"][:]
            acc = h["acc"][:]
            pe = h["potential_energy"][:]
            mass = h["mass"][:]
            ke = h["kinetic_energy"][:]
            density = h["density"][:]
            pressure = h["pressure"][:]
            smoothing_length = h["smoothing_length"][:]
            tag = h["tag"][:]

        print(f"Read state from {filename}")

        return State(mass, pos, vel, acc, pe, ke, density, pressure, smoothing_length, tag, ndim, iteration, time, name)

    @classmethod
    def empty(cls, n, ndim=3):
        time = 0.0
        iteration = 0
        ndim = ndim
        mass = np.ones(n, dtype=np.float64)
        pe = np.zeros(n, dtype=np.float64)
        ke = np.zeros(n, dtype=np.float64)
        pos = np.zeros((n, ndim), dtype=np.float64)
        vel = np.zeros((n, ndim), dtype=np.float64)
        acc = np.zeros((n, ndim), dtype=np.float64)
        density = np.zeros((n, ndim), dtype=np.float64)
        pressure = np.zeros((n, ndim), dtype=np.float64)
        smoothing_length = np.zeros((n, ndim), dtype=np.float64)

        tag = np.zeros(n, dtype=np.int32)

        return cls(mass, pos, vel, acc, pe, ke, density, pressure, smoothing_length, tag, ndim, iteration, time)

    def write(self, filename=None):
        if filename is None:
            filename = f"{self.name}_{self.iteration:05d}.h5"
        with h5py.File(filename, "w") as h:
            h.attrs["time"] = np.float64(self.time)
            h.attrs["n"] = np.uint64(len(self.mass))
            h.attrs["iteration"] = np.uint64(self.iteration)
            h.attrs["NDIM"] = np.uint64(self.pos.shape[1])
            h.attrs["name"] = self.name

            h.create_dataset("pos", data=self.pos)
            h.create_dataset("vel", data=self.vel)
            h.create_dataset("acc", data=self.acc)
            h.create_dataset("mass", data=self.mass)
            h.create_dataset("potential_energy", data=self.pe)
            h.create_dataset("kinetic_energy", data=self.ke)
            h.create_dataset("density", data=self.density)
            h.create_dataset("pressure", data=self.pressure)
            h.create_dataset("smoothing_length", data=self.smoothing_length)
            h.create_dataset("tag", data=self.tag.astype(np.int32))

        print(f"Written state to {filename}")

    def __add__(self, other):

        if self.time != other.time:
            raise ValueError("States are at different times")

        if self.iteration != other.iteration:
            raise ValueError("States are at different iterations")

        if self.ndim != other.ndim:
            raise ValueError("States are at different NDIMs")

        if self.name != other.name:
            raise ValueError("Names are not the same")

        n1 = len(self.mass)
        n2 = len(other.mass)
        n = n1 + n2
        state = State.empty(n)
        state.mass = np.hstack((self.mass, other.mass))
        state.pos = np.vstack((self.pos, other.pos))
        state.vel = np.vstack((self.vel, other.vel))
        state.acc = np.vstack((self.acc, other.acc))
        state.pe = np.hstack((self.pe, other.pe))
        state.ke = np.hstack((self.ke, other.ke))
        state.density = np.hstack((self.density, other.density))
        state.pressure = np.hstack((self.pressure, other.pressure))
        state.smoothing_length = np.hstack((self.smoothing_length, other.smoothing_length))
        state.tag = np.hstack((self.tag, other.tag))
        state.name = self.name
        state.time = self.time
        state.iteration = self.iteration

        return state
