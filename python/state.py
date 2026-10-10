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
    ndim: int
    iteration: int = 0
    time: float = 0.0

    @classmethod
    def read(cls, filename):
        with h5py.File(filename, "r") as h:

            time = h.attrs["time"]
            iteration = h.attrs["iteration"]
            ndim = h.attrs["NDIM"]

            pos = h["pos"][:]
            vel = h["vel"][:]
            acc = h["acc"][:]
            pe = h["potential_energy"][:]
            mass = h["mass"][:]
            ke = h["kinetic_energy"][:]

        print(f"Read state from {filename}")

        return State(mass, pos, vel, acc, pe, ke, ndim, iteration, time)

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

        return cls(mass, pos, vel, acc, pe, ke, ndim, iteration, time)

    def write(self, filename):
        with h5py.File(filename, "w") as h:
            h.attrs["time"] = np.float64(self.time)
            h.attrs["n"] = np.uint64(len(self.mass))
            h.attrs["iteration"] = np.uint64(self.iteration)
            h.attrs["NDIM"] = np.uint64(self.pos.shape[1])

            h.create_dataset("pos", data=self.pos)
            h.create_dataset("vel", data=self.vel)
            h.create_dataset("acc", data=self.acc)
            h.create_dataset("mass", data=self.mass)
            h.create_dataset("potential_energy", data=self.pe)
            h.create_dataset("kinetic_energy", data=self.ke)

        print(f"Written state to {filename}")
