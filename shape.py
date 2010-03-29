from scipy import conjugate, real;
from scipy.fftpack.basic import fft2, ifft2;
from scipy.fftpack.helper import fftshift, ifftshift;
from scipy.fftpack.pseudo_diffs import shift;
from polar import *;


#Converts geometry into polar fourier transform
def geom2polar_ft(real_geom, nR):
	return rect2polar(ifftshift(fft2(ifftshift(real_geom))), nR);

#Converts polar fourier transform into geometry
def polar_ft2geom(pft, xs, ys):
	return fftshift(real(ifft2(fftshift(polar2rect(pft, xs, ys)))))

#Rotates some polar geometry
def rotate_polar(pft, theta):
	R = len(pft);
	res = [];
	for k in range(R):
		res.append(shift(pft[k], theta));
	return res;

#Performs pointwise polynomial multiplication to evaluate a single point
def eval_pt(pft, x):
	return 0;
	
#Performs a convolution
def polar_conv(pX, pY):
	res = []
	for k in range(len(pX)):
		res.append(pX[k] * conjugate(pY[k]));
	return res;


