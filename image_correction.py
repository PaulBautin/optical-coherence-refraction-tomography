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
import nibabel as nib
import matplotlib.pyplot as plt
import numpy as np
from skimage.transform import PiecewiseAffineTransform, warp
import scipy.io
from skimage.transform import radon, rescale, iradon


def get_parser():
    """parser function"""
    parser = argparse.ArgumentParser(
        description="Algorithme de correction des distorsions liées aux changements d'indices de réfractions.",
        formatter_class=argparse.RawTextHelpFormatter,
        prog=os.path.basename(__file__).strip(".py")
    )

    mandatory = parser.add_argument_group("\nMANDATORY ARGUMENTS")
    mandatory.add_argument(
        "-i",
        required=True,
        default='octr_data',
        help='Path to folder that contains input images. Example: "octr_data"',
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
    """
    Creation de la carte des indices de refraction: "ri_map"  a partir d'un b_scan sur le plan 'zy'
    (z: axe de propagation laser et y: axe perpendiculaire a l'axe de rotation du capillaire)  avec indexage = 'ij'
    (i = y indice de ligne y et j = z indice de colonne)
    :param n_z: Nombre de colonne sur l'image.  Example: image.shape[1] = 576
    :param n_y: Nombre de ligne sur l'image.  Example: image.shape[0] = 512
    :param n: Dictionnaire contenant les indices de refraction des principales milieux
    (air, verre, agarose). Example: n["verre"] = 1.51
    :param centre: Tuple contenant position du centre du tube en cordonne cartésiennes 'xy'.  Example: (z_tube, y_tube) = ()
    :param rayon_int: Rayon intérieur du tube de verre
    :param rayon_ext: Rayon extérieur du tube de verre
    :return ri_map: carte des indices de refraction
    """
    # vecteur des coordonnées sur l'axe z (nombre de colonne)
    z = np.linspace(0, n_z, n_z)
    # vecteur des coordonnées sur l'axe y (nombre de ligne)
    y = np.linspace(0, n_y, n_y)
    # par defaut meshgrid utilise un indexage cartésien
    zv, yv = np.meshgrid(z, y, sparse=False, indexing='xy')
    # creation carte distance du centre
    dist_du_centre = np.sqrt((zv - centre[0])**2 + (yv - centre[1])**2)
    # masque de l'extérieur du tube capillaire
    masque_capillaire_ext = dist_du_centre <= rayon_ext
    # masque de l'intérieur du tube capillaire
    masque_capillaire_int = dist_du_centre <= rayon_int
    # initialisation de la carte d'indice de refraction indexage 'ij'
    ri_map = np.zeros([n_y, n_z])
    ri_map[masque_capillaire_ext] = n["verre"]
    ri_map[~masque_capillaire_ext] = n["air"]
    ri_map[masque_capillaire_int] = n["agarose"]
    return ri_map


def geo_tube_onclick(b_scan):
    """
    Permet de positionner les interface du capillaire sur l'image en cliquant dessus. Il faut que l'image soit fermee
    a chaque click sinon les cordonnées de l'interface ne seront pas enregistrer
    :param b_scan: slice de l'image que l'utilisateur veut corriger
    :return r_capillaire_ext: rayon extérieur du capillaire
    :return r_capillaire_int: rayon intérieur du capillaire
    :return z_tube: position du centre du tube sur l'axe z
    :return y_tube: position du centre du tube sur l'axe y
    """
    coord_list = []
    for i in range(4):
        fig, ax = plt.subplots()
        ax.imshow(b_scan, cmap="gray", origin="upper")
        cid = fig.canvas.mpl_connect('button_press_event', onclick)
        plt.show()
        coord_list.append(coords)
    # rayon extérieur capillaire
    r_capillaire_ext = np.abs(((coord_list[0][0] - coord_list[3][0]) / 2))
    # rayon a l'interieur capillaire
    r_capillaire_int = np.abs(((coord_list[1][0] - coord_list[2][0]) / 2))
    z_tube = coord_list[1][0] + r_capillaire_ext
    y_tube = (coord_list[0][1] + coord_list[1][1] + coord_list[2][1] + coord_list[3][1]) / 4
    print("le centre du tube est en (x: {}, y: {}) ".format(z_tube, y_tube))
    return r_capillaire_ext, r_capillaire_int, z_tube, y_tube

def nadaraya_watson(ri_map, pos=None, sig=1):
    """
    L'estimateur de Nadaraya-Watson est une méthode d’estimation par noyau qui permet d’estimer l’indice de
    réfraction en tout point de l’image i.e. a des position avec de indices flottant.
    pour lire plus sur l'estimateur de Nadaraya Watson: https://en.wikipedia.org/wiki/Kernel_regression
    :param ri_map: carte des indices de refraction
    :param pos: position dans la carte des indices de rerfaction. Example: (z_pos, y_pos)
    :param sig: ecart type de la fonction gaussienne (fenêtre des noyaux gaussiens)
    :return n: indice de refraction a la position pos
    :return dndz: derivee de l'indice de refraction selon z
    :return dndy: derivee de l'indicce de refraction selon y
    """
    # creation carte des indices et attribution du poids de chaque kernel en fonction de l'indice de refraction
    # et distance au points d'interet.
    n_z = ri_map.shape[1]
    n_y = ri_map.shape[0]
    z = np.linspace(0, n_z, n_z)
    y = np.linspace(0, n_y, n_y)
    zv, yv = np.meshgrid(z, y)
    one = np.ones((n_y, n_z))
    g_kernels = np.exp(-((zv - pos[0]*one) ** 2 + (yv - pos[1]*one) ** 2) / (2 * (sig*one)**2))
    norm = np.sum(g_kernels, axis=(0, 1))
    n_A = np.sum((ri_map * g_kernels), axis=(0, 1))
    # indice de refraction en la position d'interet
    n = n_A / norm

    # calcul analytique des derivees
    # calcul analytique de la derive du numérateur en z et y a la position d'interet
    dn_Adz = ri_map * g_kernels * (zv - pos[0]*one) / sig ** 2
    dn_Ady = ri_map * g_kernels * (yv - pos[1]*one) / sig ** 2
    dn_Adz = np.sum(dn_Adz, axis=(0, 1))
    dn_Ady = np.sum(dn_Ady, axis=(0, 1))
    # calcul analytique de la derive du denominateur en z et y a la position d'interet
    dnormdz = g_kernels * (zv - pos[0]*one) / sig ** 2
    dnormdy = g_kernels * (yv - pos[1]*one) / sig ** 2
    dnormdz = np.sum(dnormdz, axis=(0, 1))
    dnormdy = np.sum(dnormdy, axis=(0, 1))
    # formule de derive d'une fraction
    dndz = (norm * dn_Adz - n_A * dnormdz) / norm ** 2
    dndy = (norm * dn_Ady - n_A * dnormdy) / norm ** 2
    return n, dndz, dndy


def eq_rayon(z, y, f, ri_map):
    """
    equation differentielle resolue avec la methode RK4
    :param z: position du rayon sur l'axe z (z: axe propagation rayon)
    :param y: position du rayon sur l'axe y (z: axe perpendiculaire a la propagation rayon)
    :param f: paramètre servant a résoudre RK4 d'une equation diff du second degree
    :param ri_map: carte des indices de refraction
    :return: dfdz = d^2y/dz^2
    """
    # equation differentielle a rentrer dans RK4
    n, dndz, dndy = nadaraya_watson(ri_map, pos=[z, y])
    dfdz = 1. / n * (dndy - dndz * f) * (1 + f ** 2)
    return dfdz


def rk4(ri_map):
    """
    Fonction permettant la resolution numerique de la propagation des rayons optiques dans un milieux avec
    des indices de refractions non-homogene. La resolution numerique est basee sur la methode de Runge-Kutta
    d'ordre 4. Les condition au frontiere sont: dydz(0)=0, z=0 . Chaque pas de rayons est divise par l'indice
    de refraction du milieux.
    :param ri_map: carte des indices de refraction
    :return src: liste des indices positionnelles des rayons sans changement d'indice de refraction
    (Coordonnees de l'image source)
    :return dst: liste des indices positionelles des rayons dans un milieux d'indice de refraction non homogeme
    (Coordonnees de l'image destination)
    """
    # programmation RK4
    n_y, n_z = ri_map.shape
    src_rows = []
    src_cols = []
    dst_rows = []
    dst_cols = []
    for yi in np.linspace(1, n_y-1, 7):
        fi = 0
        zi = 0
        t0 = 10
        t = 10
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
    """
    Fonction permettant de positionner les interfaces du capillaire avec un clic de la souris
    :param event:
    :return:
    """
    global coords
    coords = [event.xdata, event.ydata]
    print('position interface, xdata={}, ydata={}'.format(event.xdata, event.ydata))


def radon_f(b_scans):
    # Variable comprenant les angles des B-scans (pris dans l'article original 1 B-scan par 6 deg sur 360 deg)
    theta = np.linspace(0, 360, len(b_scans))

    # Afficher l'image pour l'angle 0
    plt.title("B-scan pour angle 0")
    plt.imshow(b_scans[0])
    plt.show()

    # Évaluer la projection pour l'angle 0 pour initialiser la variable projection
    projections = radon(b_scans[0], [0])
    # Empiler les projections pour chaque angle pour créer le sinogramme
    for i in range(len(b_scans)-1):
        projections = np.hstack((projections, radon(b_scans[i+1], [0])))
    print("les dimensions du sinogramme sont:  ".format(projections.shape))

    # Afficher le sinogramme
    plt.title("Radon transform\n(Sinogram)")
    plt.xlabel("Projection axis")
    plt.ylabel("Intensity")
    plt.imshow(projections)
    plt.show()

    # Transformee de radon inverse avec un filtre de hanning
    reconstruction = iradon(projections, theta,  filter_name='ramp', interpolation='cubic')
    # Afficher le resultat de la transformeee de radon inverse
    plt.title("Reconstruction\nfrom sinogram")
    plt.imshow(reconstruction, cmap=plt.cm.Greys_r)
    plt.show()


def main():
    """
    Fonction main,
    """
    # get parser elements
    parser = get_parser()
    arguments = parser.parse_args()
    path_data = os.path.abspath(os.path.expanduser(arguments.i))
    path_output = os.path.abspath(arguments.o)

    # get image (z: axe de propagation laser, y: axe de rotation capillaire, x: axe perpendiculaire
    # a la rotation capillaire)
    filename, file_extension = os.path.splitext("insect_leg_data.mat")
    path_image = os.path.join(path_data, filename + file_extension)
    print("L'image est prise dans le fichier {}".format(path_image))
    # Les images sur le Github du papier original sont en format Matlab (extension: .mat), il faut donc les transformer
    if file_extension == ".mat":
        mat = scipy.io.loadmat(path_image)
        print(mat.keys())
        nz, ny = mat['Bscan_dims'][0][0], mat['Bscan_dims'][0][1]
        print(len(mat["xzcoords"]))
        print(nz)
        print(ny)
        b_scans = mat["Bscans"]
        b_scan = b_scans[0]
    else:
        image = nib.load(path_image)
        nx, ny, nz = image.shape
        print("Le format de l'image est (x:{}, y:{}, z:{})".format(nx, ny, nz))
        # extraire un b_scan avec plan 'zy'
        data = image.get_fdata()
        x_tube = 310
        b_scan = data[:, x_tube, :]

    radon_f(b_scans)

    # creation de la carte des indices de refraction a partir de la geometrie du tube
    # geometrie tube fabricant:  Longueur : 75 mm. Diamètre intérieur : 1,1 à 1,2 mm. Paroi : 0,2 mm ± 0,02 mm
    # argument.l sert a rentrer la position des interfaces avec un click de la souris
    if arguments.l:
        r_capillaire_ext, r_capillaire_int, z_tube, y_tube = geo_tube_onclick(b_scan)
    # sinon les coordonnees seront rentree manuellement ci-dessous
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
    # creation de la carte des indices de refraction et affichage
    ri_map_init = ri_map(nz, ny, n, (z_tube, y_tube), r_capillaire_int, r_capillaire_ext)
    plt.imshow(ri_map_init, cmap="gray")
    plt.show()

    # tracer les rayons avec la fonction rk4
    src, dst = rk4(ri_map_init)
    plt.imshow(ri_map_init, cmap="gray")
    plt.show()
    # np.save permet de sauver la liste des coordonnes de chaque rayon
    #np.save('source_coord', src)
    #np.save('destination_coord', dst)
    #src = np.load('source_coord.npy')
    #dst = np.load('destination_coord.npy')

    # "Control points are used to define the mapping. The transform is based on a Delaunay triangulation
    # of the points to form a mesh. Each triangle is used to find a local affine transform"
    # https://scikit-image.org/docs/dev/api/skimage.transform.html?highlight=transform#skimage.transform.PiecewiseAffineTransform
    tform = PiecewiseAffineTransform()
    tform.estimate(dst, src)
    # warp b_scan to new coordinates
    out_rows = b_scan.shape[0]
    out_cols = b_scan.shape[1]
    out = warp(b_scan, tform, output_shape=(out_rows, out_cols))
    # plot reconstructed b_scan
    fig, ax = plt.subplots()
    ax.imshow(out)
    #ax.plot(tform.inverse(src)[:, 0], tform.inverse(src)[:, 1], '.b')
    ax.axis((0, out_cols, out_rows, 0))
    if arguments.o is not None:
        filename_out = os.path.splitext(os.path.basename(filename))[0] + '_recon.png'
        plt.savefig(os.path.join(path_output, filename_out))
    plt.show()


if __name__ == "__main__":
    main()