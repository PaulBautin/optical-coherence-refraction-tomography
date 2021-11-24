
import scipy.io
from scipy import interpolate
from scipy import signal
import nibabel as nib
import matplotlib.pyplot as plt
import numpy as np
from skimage.transform import PiecewiseAffineTransform, warp
from skimage import data
import math

# parametres
# geometrie du tube en mm
# fabricant:  Longueur : 75 mm. Diamètre intérieur : 1,1 à 1,2 mm. Paroi : 0,2 mm ± 0,02 mm
d_int_tube = 1.1
d_ext_tube = 1.3
# pixel size in mm
pix_dim = 1.1 / 160
r_capillaire_ext = ((362 - 115)/2)
r_capillaire_int = ((319 - 159)/2)
z_tube = 115 + r_capillaire_ext
x_tube = 260
y_tube = 267

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
b_scan = data[x_tube, :, :]
plt.imshow(b_scan, cmap="gray", origin="lower")
plt.show()

# carte des indices de refraction (ri_map)
def ri_map(n_x, n_y, n, centre=None, rayon_int=None, rayon_ext=None):
    x = np.linspace(0, n_x, n_x)
    y = np.linspace(0, n_y, n_y)
    xv, yv = np.meshgrid(x, y)
    dist_du_centre = np.sqrt((xv - centre[0])**2 + (yv - centre[1])**2)
    masque_capillaire_ext = dist_du_centre <= rayon_ext
    masque_capillaire_int = dist_du_centre <= rayon_int
    ri_map = np.zeros([n_y, n_x])
    ri_map[masque_capillaire_ext] = n["verre"]
    ri_map[~masque_capillaire_ext] = n["air"]
    ri_map[masque_capillaire_int] = n["agarose"]
    return ri_map

def nadaraya_watson(ri_map, pos=None, sig=1):
    n_z = ri_map.shape[1]
    n_y = ri_map.shape[0]
    z = np.linspace(0, n_z, n_z)
    y = np.linspace(0, n_y, n_y)
    zv, yv = np.meshgrid(z, y)
    one = np.ones((n_y, n_z))
    g_kernels = np.exp(-((zv - pos[0]*one) ** 2 + (yv - pos[1]*one) ** 2) / (2 * (sig*one)**2))
    norm = np.sum(g_kernels, axis=(0, 1))
    n_A = np.sum((ri_map * g_kernels), axis=(0, 1))

    dn_Adz = ri_map * g_kernels * (zv - pos[0]*one) / sig ** 2
    dn_Ady = ri_map * g_kernels * (yv - pos[1]*one) / sig ** 2
    dn_Adz = np.sum(dn_Adz, axis=(0, 1))
    dn_Ady = np.sum(dn_Ady, axis=(0, 1))


    dnormdz = g_kernels * (zv - pos[0]*one) / sig ** 2
    dnormdy = g_kernels * (yv - pos[1]*one) / sig ** 2
    dnormdz = np.sum(dnormdz, axis=(0, 1))
    dnormdy = np.sum(dnormdy, axis=(0, 1))

    n = n_A / norm
    dndz = (norm * dn_Adz - n_A * dnormdz) / norm ** 2
    dndy = (norm * dn_Ady - n_A * dnormdy) / norm ** 2
    return n, dndz, dndy

def nadaraya_watson_index(value, object_map, pos=None, sig=0.5):
    n_z = object_map.shape[1]
    n_y = object_map.shape[0]
    z = np.linspace(0, n_z, n_z)
    y = np.linspace(0, n_y, n_y)
    zv, yv = np.meshgrid(z, y)
    one = np.ones((n_y, n_z))
    object = one * value
    g_kernels = np.exp(-((zv - pos[0]*one) ** 2 + (yv - pos[1]*one) ** 2) / (2 * (sig*one)**2))
    object = object * g_kernels
    return object

def eq_rayon(z, y, f, ri_map):
    # equation differentiel a rentrer dans RK4
    #print("position {} et {}".format(z, y))
    n, dndz, dndy = nadaraya_watson(ri_map_init, pos=[z, y])
    dfdz = 1. / n * (dndy - dndz * f) * (1 + f ** 2)
    return dfdz


def rk4(ri_map, b_scan, start):
    # programmation RK4 dydz(0)=0 y(0)=1
    n_y, n_z = ri_map.shape
    object_map = np.zeros((n_y, n_z))
    for yi in range(start[1], start[1]+401, 25):
        start = [0, yi]
        fi = 0
        n, _, _ = nadaraya_watson(ri_map, pos=start)
        t = 1
        step = t / n
        object_map[yi, 0] = b_scan[yi, 0]
        zi = 0
        z_n = []
        y_n = []
        while (0 <= yi < n_y) & (0 <= zi < n_z) & (0 <= t < n_z-1):
            print("le code avance yi = {} zi = {}".format(yi, zi))
            l0 = step * eq_rayon(zi, yi, fi, ri_map)
            k0 = fi

            z1 = zi + step / 2
            y1 = yi + k0 / 2
            f1 = fi + l0 / 2
            l1 = step * eq_rayon(z1, y1, f1, ri_map)
            k1 = step * f1

            z2 = zi + step / 2
            y2 = yi + k1 / 2
            f2 = fi + l1 / 2
            l2 = step * eq_rayon(z2, y2, f2, ri_map)
            k2 = step * f2

            z3 = zi + step / 2
            y3 = yi + k2 / 2
            f3 = fi + l2 / 2
            l3 = step * eq_rayon(z3, y3, f3, ri_map)
            k3 = step * f3

            zi = zi + step
            yi = yi + (k0 + 2 * k1 + 2 * k2 + k3) / 6
            fi = fi + (l0 + 2 * l1 + 2 * l2 + l3) / 6
            t = t + 1
            z_n.append(zi)
            y_n.append(yi)
            object_map = object_map + nadaraya_watson_index(b_scan[start[1], t], object_map, pos=[zi, yi], sig=1)
        # print("iteration = {}   ,  y = {}".format(zi, yi))
        plt.plot(z_n, y_n, "r")
    print("sortie de la boucle {}".format(yi))
    return object_map

ri_map_init = ri_map(b_scan.shape[1], b_scan.shape[0], n, (z_tube, y_tube), r_capillaire_int, r_capillaire_ext)
plt.imshow(ri_map_init, cmap="gray")
plt.show()
plt.imshow(ri_map_init, cmap="gray")
object_map = rk4(ri_map_init, b_scan, start=[0, y_tube-200])
plt.show()
plt.imshow(object_map, cmap="gray")
plt.show()


img = b_scan
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


fig, ax = plt.subplots()
ax.imshow(out)
ax.plot(tform.inverse(src)[:, 0], tform.inverse(src)[:, 1], '.b')
ax.axis((0, out_cols, out_rows, 0))
plt.show()
