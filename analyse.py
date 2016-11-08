# note: depends on a scipy and matplotlib
#       easiest way to get them, is to install Anaconda, https://www.continuum.io/anaconda-overview

import sys
import scipy
import scipy.misc
import numpy as np
from matplotlib import pyplot as plt

if ( len( sys.argv ) ) < 2:
	sys.exit("usage: analyse.py <filename.bmp>")

#########
# init
img = scipy.misc.imread( sys.argv[1], mode='L' )

###########
# calc
f = np.fft.fft2(img)

# print( f )

fshift = np.fft.fftshift(f)
# magnitude_spectrum = 20*np.log(np.abs(fshift))
magnitude_spectrum = np.log(np.abs(fshift))
phase_spectrum = np.angle(fshift)

# TODO: average radial samples and plot as graph (either sample radially, or do transform and average rows)
# TODO: some way to plot initial thresholded "animation"? (consequtive points)

###########
# plotting

#TODO: use extent to force pixel-perfect?

interp = 'bicubic' #nearest

# plt.clf()

plt.subplot(121)
plt.imshow(img, cmap = 'gray', interpolation=interp)
plt.title('Input Image:\n' + sys.argv[1] ), plt.xticks([]), plt.yticks([])

plt.subplot(122)
# plt.imshow(magnitude_spectrum, interpolation=interp)
plt.imshow(magnitude_spectrum, cmap = 'gray', interpolation=interp)
# plt.imshow(magnitude_spectrum, cmap = 'gray', interpolation=interp, vmin=8, vmax=12)
plt.title('Magnitude Spectrum'), plt.xticks([]), plt.yticks([])

# plt.subplot(133)
# plt.imshow(phase_spectrum, cmap='gray', interpolation=interp)
# plt.title('Phase Spectrum'), plt.xticks([]), plt.yticks([])

plt.savefig( 'figure_' + sys.argv[1] + '.png' )
# plt.show()
