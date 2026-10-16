"""
Example script to demonstrate creating a galaxy and visualising it
"""
import numpy as np
import matplotlib.pyplot as plt

from tools.profiles import ExponentialProfile, HernquistProfile, IsothermalProfile
from tools.galaxy_model import Isotropic3D, Disk3D, Galaxy
from tools.state import State
from tools.visualise import make_image

# Create individual components of galaxy (here a disk and a bulge)
disk = Disk3D(10000, ExponentialProfile(1, 0.3), ExponentialProfile(0.2, 0.02), 1, tag=1)
bulge = Isotropic3D(1000, HernquistProfile(0.2, 0.1), 0.2, tag=2)
halo = Isotropic3D(2000, IsothermalProfile(1, 0.05), 4, tag=-1)

# Create galaxy from these components
gal = Galaxy.from_components([disk, bulge, halo])
# Rotate
gal.rotate(30, 60)
gal.set_experiment_name("example")
gal.write()

# N-Body code can now be run using this file as an input

# Read in state and visualise it
state = State.read("example_00000.h5")
x, y, z = state.pos.T
img = make_image(x, z, l0=2.5)
plt.imshow(np.arcsinh(img), cmap="bone")
plt.show()
