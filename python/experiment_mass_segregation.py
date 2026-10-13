from profiles import PlummerProfile
from galaxy_model import GalaxyComponent, Isotropic3D, Galaxy
import matplotlib.pyplot as plt

m1 = Isotropic3D(1000, PlummerProfile(1, 0.2), mass=1.0)
m2 = Isotropic3D(1000, PlummerProfile(1, 0.2), mass=2.0)
m5 = Isotropic3D(1000, PlummerProfile(1, 0.2), mass=5.0)
m10 = Isotropic3D(1000, PlummerProfile(1, 0.2), mass=10.0)

gal = Galaxy.from_components([m1, m2, m5, m10])
gal.normalise_masses_for_time()

x, y, z = gal.state.pos.T
u, v, w = gal.state.vel.T

plt.quiver(x, y, u, v)
plt.show()

gal.write("massseg_00000.h5")
