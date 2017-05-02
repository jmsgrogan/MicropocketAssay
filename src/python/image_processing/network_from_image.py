from skimage.data import camera
from skimage.filters import frangi, hessian
from skimage.morphology import skeletonize
from skimage.color import rgb2grey
from skimage.filters import threshold_minimum
from skimage import io
import scipy.stats
import numpy as np

import matplotlib.pyplot as plt


current_image = io.imread("/home/grogan/test2.tif")
#cropped = image[crop_top:crop_bottom]
fig, ax = plt.subplots()
#ax.imshow(frangi(current_image))

# Crop to ROI
#crop = current_image[260:380, 180:600]
crop = current_image[300:480, 250:650]

# Get luminance from RGB
lum = rgb2grey(crop)

# Remove high luminance artifacts
print np.max(lum)
thres = np.clip(lum, 0, 0.47)

# Invert so vessels are light and rescale
inv = thres

edges = frangi(inv, scale_range = (0.0, 4.0), scale_step = 0.5, black_ridges=True)

print np.max(edges)
binary_min = scipy.stats.threshold(edges, threshmin=0, threshmax=3.e-7, newval=1)

skeleton = skeletonize(binary_min)
#ax.imshow(thres)
ax.imshow(edges)
ax.imshow(binary_min)
ax.imshow(skeleton)

#io.imsave("/home/grogan/frangi2.tif", edges)

plt.show()