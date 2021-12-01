# !/usr/bin/env python
# -*- coding: utf-8
#########################################################################################
#
# Corriger les distortions sur l'image du au changement d'indice de refraction
#
# ---------------------------------------------------------------------------------------
# Authors: groupe oct
#
#########################################################################################


import os
import argparse
import scipy.io
from scipy import interpolate
from scipy import signal
import nibabel as nib
import matplotlib.pyplot as plt
import numpy as np
from skimage.transform import PiecewiseAffineTransform, warp
from skimage import data
import math


def get_parser():
    """parser function"""
    parser = argparse.ArgumentParser(
        description="Correction de l'image par rapport au changment d'indice de refraction",
        formatter_class=argparse.RawTextHelpFormatter,
        prog=os.path.basename(__file__).strip(".py")
    )

    mandatory = parser.add_argument_group("\nMANDATORY ARGUMENTS")
    mandatory.add_argument(
        "-i",
        required=True,
        default='octr_data',
        help='Path to folder that contains input images (e.g. "octr_data")',
    )
    optional = parser.add_argument_group("\nOPTIONAL ARGUMENTS")
    optional.add_argument(
        '-fig',
        help='Generate figures',
        action='store_true'
    )
    optional.add_argument(
        '-o',
        help='Path where figures will be saved. By default, they will be saved in the current directory.',
        default=""
    )
    optional.add_argument(
        '-l',
        help='manually find borders of capillary',
        action='store_true'
    )
    return parser


def ri_map(n_z, n_y, n, centre=None, rayon_int=None, rayon_ext=None):
    """ Creation de la carte des indices de refraction: "ri_map"  a partir d'un b_scan sur le plan 'zy'
    (z: axe de propagation laser et y: axe perpenticulaire a l'axe de rotation du capilllaire)  avec indexage = 'ij'
    (i = y indice de ligne y et j = z indice de colonne)
    :param n_z: Nombre de colonne sur l'image.  Example: image.shape[1] = 576
    :param n_y: Nombre de ligne sur l'image.  Example: image.shape[0] = 512
    :param n: Dictionnaire contennant les indices de refraction des principales milieux
    (air, verre, agarose). Example: n["verre"] = 1.51
    :param centre: Tuple contennant position du centre du tube cordonne cartesienne 'xy'.  Example: (z_tube, y_tube) = ()
    :param rayon_int: Rayon interieur du tube de verre
    :param rayon_ext: Rayon exterieur du tube de verre
    :return ri_map: carte des indices de refraction
    """
    # vecteur des coordonnees sur l'axe z (nombre de colonne)
    z = np.linspace(0, n_z, n_z)
    # vecteur des coordonnees sur l'axe y (nombre de ligne)
    y = np.linspace(0, n_y, n_y)
    # par defaut meshgrid utilise un indexage cartesien
    zv, yv = np.meshgrid(z, y, sparse=False, indexing='xy')
    # creation carte distance du centre
    dist_du_centre = np.sqrt((zv - centre[0])**2 + (yv - centre[1])**2)
    # masque de l'exterieur du tube capillaire
    masque_capillaire_ext = dist_du_centre <= rayon_ext
    # masque de l'interieur du tube capillaire
    masque_capillaire_int = dist_du_centre <= rayon_int
    # initialisation carte de d'indice de refraction indexage 'ij'
    ri_map = np.zeros([n_y, n_z])
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
    n, dndz, dndy = nadaraya_watson(ri_map, pos=[z, y])
    dfdz = 1. / n * (dndy - dndz * f) * (1 + f ** 2)
    return dfdz


def rk4(ri_map):
    # programmation RK4 dydz(0)=0 y(0)=1
    n_y, n_z = ri_map.shape
    src_rows = []
    src_cols = []
    dst_rows = []
    dst_cols = []
    #object_map = np.zeros((n_y, n_z))
    for yi in np.linspace(1, n_y-1, 10):
        fi = 0
        zi = 0
        t0 = 5
        t = 5
        t_i = [0]
        y_i = [yi]
        z_n = [0]
        y_n = [yi]
        while (0 <= yi < n_y) & (0 <= zi < n_z) & (0 <= t < n_z-1):
            n, _, _ = nadaraya_watson(ri_map, pos=[zi, yi])
            step = t0 / n
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
            t += t0
            t_i.append(t)
            y_i.append(y_n[0])
            z_n.append(zi)
            y_n.append(yi)
        # print("iteration = {}   ,  y = {}".format(zi, yi))
        src_rows = src_rows + y_i
        src_cols = src_cols + t_i
        dst_rows = dst_rows + y_n
        dst_cols = dst_cols + z_n
        plt.plot(z_n, y_n, "r")
        plt.plot(t_i, y_i, "b")
    src = np.vstack([src_cols, src_rows]).T
    dst = np.vstack([dst_cols, dst_rows]).T
    return src, dst


def onclick(event):
    global coords
    coords = [event.xdata, event.ydata]
    print('position interface, xdata={}, ydata={}'.format(event.xdata, event.ydata))


def main():
    """
    main function, gather stats and call plots
    """
    # get parser elements
    parser = get_parser()
    arguments = parser.parse_args()
    path_data = os.path.abspath(os.path.expanduser(arguments.i))
    print("L'image est prise dans le fichier {}".format(path_data))
    path_output = os.path.abspath(arguments.o)

    # get image (z: axe de propagation laser, y: axe de rotation capillaire, x: axe perpendiculaire a la rotation capillaire)
    filename = "data_test_oct_0degree.nii"
    image = nib.load(os.path.join(path_data, filename))
    nx, ny, nz = image.shape
    print("Le format de l'image est (x:{}, y:{}, z:{})".format(nx, ny, nz))

    # extraire un b_scan avec plan 'zy'
    data = image.get_fdata()
    x_tube = 260
    b_scan = data[x_tube, :, :]

    # geometrie du tube
    # fabricant:  Longueur : 75 mm. Diamètre intérieur : 1,1 à 1,2 mm. Paroi : 0,2 mm ± 0,02 mm
    if arguments.l:
        coord_list = []
        for i in range(4):
            fig, ax = plt.subplots()
            ax.imshow(b_scan, cmap="gray", origin="upper")
            cid = fig.canvas.mpl_connect('button_press_event', onclick)
            plt.show()
            coord_list.append(coords)
        # rayon exterieur capillaire
        r_capillaire_ext = np.abs(((coord_list[0][0] - coord_list[3][0]) / 2))
        # rayon a l'interieur capillaire
        r_capillaire_int = np.abs(((coord_list[1][0] - coord_list[2][0]) / 2))
        z_tube = coord_list[1][0] + r_capillaire_ext
        y_tube = (coord_list[0][1] + coord_list[1][1] + coord_list[2][1] + coord_list[3][1]) / 4
        print("le centre du tube est en (x: {}, y: {}) ".format(z_tube, y_tube))
    else:
        r_capillaire_ext = ((362 - 115)/2)
        r_capillaire_int = ((319 - 159)/2)
        z_tube = 115 + r_capillaire_ext
        y_tube = 267

    # propriete optique
    n = {
        "air": 1.000293,
        "verre": 1.51,
        "agarose": 1.34,
    }

    # plot le la carte des indices de refraction
    ri_map_init = ri_map(nz, ny, n, (z_tube, y_tube), r_capillaire_int, r_capillaire_ext)
    plt.imshow(ri_map_init, cmap="gray")
    plt.show()

    # tracer les rayons avec la fonction rk4
    #src, dst = rk4(ri_map_init)
    #plt.imshow(ri_map_init, cmap="gray")
    #plt.show()

    #np.save('source_coord', src)
    #np.save('destination_coord', dst)

    src = np.load('source_coord.npy')
    dst = np.load('destination_coord.npy')

    # plot l'image corrigee
    # plt.imshow(object_map, cmap="gray")
    # plt.show()
    #print(np.max(src - dst))
    #print(dst.shape)

    tform = PiecewiseAffineTransform()
    tform.estimate(dst, src)

    out_rows = b_scan.shape[0]
    out_cols = b_scan.shape[1]
    out = warp(b_scan, tform, output_shape=(out_rows, out_cols))

    fig, ax = plt.subplots()
    ax.imshow(out)
    #ax.plot(tform.inverse(src)[:, 0], tform.inverse(src)[:, 1], '.b')
    ax.axis((0, out_cols, out_rows, 0))
    plt.show()


if __name__ == "__main__":
    main()