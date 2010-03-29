from scipy import *;
from scipy.fftpack.fft
from polar import *;


#Converts geometry into polar fourier transform
def geom2polar_ft(real_geom, nR):
	return rect2polar(ifftshift(rfft2(ifftshift(real_geom))), nR);

#Converts polar fourier transform into geometry
def polar_ft2geom(pft, xs, ys):
	return fftshift(irfft2(fftshift(polar2rect(polar_ft, xs, ys))))

#Performs a circular rotation
def circ_rotate(circ, theta):
	return circ;

#Rotates some polar geometry
def rotate_polar(pft, theta):
	R = len(pft);
	res = [];
	for k in range(R):
		res.append(circ_rotate(pft[k], theta));
	return R;

