from scipy import conjugate, real, array;
from scipy.linalg.basic import norm;
from scipy.fftpack.basic import fft2, ifft2;
from scipy.fftpack.helper import fftshift, ifftshift;
from scipy.fftpack.pseudo_diffs import shift;
from polar import *;


#Converts geometry into polar fourier transform
def geom2polar_ft(real_geom, nR):
	return rect2polar(ifftshift(fft2(ifftshift(__pad(real_geom)))), nR);

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
	
#Pads geometry so that it doesn't have overlap issues
def __pad(f):
	ns = int(round(norm(array(f.shape))));
	res = zeros((ns,ns));
	c0 = array((ns / 2. - array(f.shape)/2.).round(), dtype('int'));
	c1 = c0 + array(f.shape, dtype('int'));
	res[c0[0]:c1[0], c0[1]:c1[1]] = f;
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


