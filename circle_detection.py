import numpy as np
import matplotlib.pyplot as plt
from skimage.transform import hough_circle, hough_circle_peaks
from skimage.feature import canny
import nibabel as nib

# load image and find center b_scan
filename = "data_test_oct_0degree.nii"
image = nib.load(filename)
data = image.get_fdata()
x_tube = 260
b_scan = data[x_tube, :, :]
print("la detection de la frontiere circulaire peu prendre du temps...")

# Load picture and detect edges
edges = canny(b_scan, sigma=2) #, low_threshold=50, high_threshold=100)
plt.imshow(edges, cmap="gray", origin="lower")
plt.show()

# Detect one radii
hough_radii = np.arange(100, 140, 1)
hough_res = hough_circle(edges, hough_radii)

# Select the most prominent circle and detect center
accums, cx, cy, radii = hough_circle_peaks(hough_res, hough_radii,
                                           total_num_peaks=1)
print("Centre du cercle sur l'axe x: {} et sur l'axe y: {}".format(cx, cy))


# Draw circle
plt.imshow(b_scan, cmap="gray", origin="lower")
circle1 = plt.Circle((cx, cy), radii, color='b', fill=False)
plt.gca().add_patch(circle1)
plt.show()