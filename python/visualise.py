import numpy as np


def make_image(x, y, masses, l0=None, xc=None, yc=None, size=1000, psize=10):
    """Make an image of the particles, where each particle is represented by a Gaussian of width psize"""

    if l0 is None:
        l0 = max(x.max() - x.min(), y.max() - y.min())

    if psize is None:
        psize = size * 0.01

    if xc is None or yc is None:
        xc = (x * masses).sum() / masses.sum()
        yc = (y * masses).sum() / masses.sum()

    # construct PSF of particles
    stampsize = int(psize) * 5
    stampsize += stampsize % 2
    kx = ky = np.fft.fftfreq(stampsize)

    fker = np.outer(np.exp(-(ky**2) * 2 * psize**2), np.exp(-(kx**2) * 2 * psize**2))

    stamp = np.fft.ifft2(fker).real
    stamp = np.fft.fftshift(stamp)
    stamp /= stamp.max()

    img = np.zeros((size + 2 * stampsize, size + 2 * stampsize))

    for i in range(len(x)):
        xp = (x[i] - xc) / l0 * size + size / 2
        yp = (y[i] - yc) / l0 * size + size / 2

        # particle outside of the image
        if xp > size or xp < 0 or yp > size or yp < 0:
            continue

        xp += stampsize
        yp += stampsize

        xp = int(xp)
        yp = int(yp)

        img[yp - stampsize // 2 : yp + stampsize // 2, xp - stampsize // 2 : xp + stampsize // 2] += stamp * masses[i]

    img = img[stampsize:-stampsize, stampsize:-stampsize]

    return img
