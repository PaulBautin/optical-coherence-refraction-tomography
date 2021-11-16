
import scipy.io
from scipy import interpolate
from scipy import signal
import nibabel as nib
import matplotlib.pyplot as plt
import numpy as np
from skimage.transform import PiecewiseAffineTransform, warp
from skimage import data

# parametres
# geometrie du tube en mm
# fabricant:  Longueur : 75 mm. Diamètre intérieur : 1,1 à 1,2 mm. Paroi : 0,2 mm ± 0,02 mm
d_int_tube = 1.1
d_ext_tube = 1.3
# pixel size in mm
pix_dim = 1.1 / 160
r_capillaire_ext = ((362 - 115)/2)*pix_dim
r_capillaire_int = ((319 - 159)/2)*pix_dim
z_tube_ext = (115 + (362 - 115)/2)*pix_dim
z_tube_int = (159 + (319 - 159)/2)*pix_dim
z_tube = (z_tube_int + z_tube_ext)/2
x_tube = 260*pix_dim
y_tube = 267*pix_dim

# propriete optique
n = {
    "air": 1.000293,
    "verre": 1.51,
    "agarose": 1.34,
}



#propriete image en pixels
filename = "data_test_oct_0degree.nii"
image = nib.load(filename)
data = image.get_fdata()
b_scan = data[round(x_tube/pix_dim), :, :]
plt.imshow(b_scan, cmap="gray", origin="lower")
plt.show()

# carte des indices de refraction (ri_map)
def ri_map(n_x, n_y, pix_dim, n, centre=None, rayon_int=None, rayon_ext=None):
    x = np.linspace(0, n_x * pix_dim, n_x)
    y = np.linspace(0, n_y * pix_dim, n_y)
    xv, yv = np.meshgrid(x, y)
    dist_du_centre = np.sqrt((xv - centre[0])**2 + (yv - centre[1])**2)
    masque_capillaire_ext = dist_du_centre <= rayon_ext
    masque_capillaire_int = dist_du_centre <= rayon_int
    ri_map = np.zeros([n_y, n_x])
    ri_map[masque_capillaire_ext] = n["verre"]
    ri_map[~masque_capillaire_ext] = n["air"]
    ri_map[masque_capillaire_int] = n["agarose"]
    return ri_map


def eq_rayon(z, y, f, ri_map):
    # equation differentiel a rentrer dans RK4
    dndz = np.gradient(ri_map)[0]
    dndy = np.gradient(ri_map)[1]
    dfdz = 1. / ri_map[y, z] * (dndy[y, z] - dndz[y, z] * f) * (1 + f ** 2)
    return dfdz


def rk4(ri_map):
    # programmation RK4 dydz(0)=0 y(0)=1
    n_y, n_z = ri_map_init.shape
    h = 2
    zi = 0
    fi = 0
    for yi in range(260, 511, 260):
        z_n = []
        y_n = []
        zi = 0
        fi = 0
        while (zi <= n_z-1) & (yi <= n_y-1):
            l0 = h * eq_rayon(zi, yi, fi, ri_map)
            k0 = fi

            z1 = int(zi + h / 2)
            y1 = int(yi + k0 / 2)
            f1 = fi + l0 / 2
            l1 = h * eq_rayon(z1, y1, f1, ri_map)
            k1 = h * f1

            z2 = int(zi + h / 2)
            y2 = int(yi + k1 / 2)
            f2 = fi + l1 / 2
            l2 = h * eq_rayon(z2, y2, f2, ri_map)
            k2 = h * f2

            z3 = int(zi + h / 2)
            y3 = int(yi + k2 / 2)
            f3 = fi + l2 / 2
            l3 = h * eq_rayon(z3, y3, f3, ri_map)
            k3 = h * f3

            zi = zi + h
            yi = round(yi + (k0 + 2 * k1 + 2 * k2 + k3) / 6)
            fi = fi + (l0 + 2 * l1 + 2 * l2 + l3) / 6
            z_n.append(zi)
            y_n.append(yi)
        plt.plot(z_n, y_n)
        print(z_n)
        print("iteration = {}   ,  y = {}".format(zi, yi))
    return z_n, y_n

ri_map_init = ri_map(b_scan.shape[1], b_scan.shape[0], pix_dim, n, (z_tube, y_tube), r_capillaire_int, r_capillaire_ext)
z_n, y_n = rk4(ri_map_init)
y_n = np.array(y_n)
plt.imshow(ri_map_init, cmap="gray")
plt.show()


img = b_scan
print(img.shape)
rows, cols = img.shape[0], img.shape[1]

src_cols = np.linspace(0, cols, 20)
src_rows = np.linspace(0, rows, 20)
src_rows, src_cols = np.meshgrid(src_rows, src_cols)
src = np.dstack([src_cols.flat, src_rows.flat])[0]

# add sinusoidal oscillation to row coordinates
y_n = signal.resample(y_n, src.shape[0])
dst_rows = src[:, 1] - (y_n - y_n[0])
dst_cols = src[:, 0]
dst_rows *= 1.5
dst_rows -= 1.5 * 50
dst = np.vstack([dst_cols, dst_rows]).T


tform = PiecewiseAffineTransform()
tform.estimate(src, dst)

out_rows = img.shape[0] - 1.5 * 50
out_cols = cols
out = warp(img, tform, output_shape=(out_rows, out_cols))
print(out.shape)

fig, ax = plt.subplots()
ax.imshow(out)
ax.plot(tform.inverse(src)[:, 0], tform.inverse(src)[:, 1], '.b')
ax.axis((0, out_cols, out_rows, 0))
plt.show()
