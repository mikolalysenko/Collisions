from scipy import *
from numpy import *
from scipy.misc import *
from scipy.fftpack import *
from scipy.signal import *
from polar import *
from se2 import *

def se2_conv(f, g, R):
	r = np.zeros((f.shape[0], f.shape[1], R))
	for i in range(R):
		r[:,:,i] = real(fftconvolve(imrotate(f, i / float(R) * 360.), g, 'same'))
	return r

def diff_polar(f, dtheta, dr):
	return real(ipfft( pft_polar_deriv(pfft(f), dtheta, dr) ))

def diff_cartesian(f, dx, dy):
	ft = fftshift(fft2(f))
	c = f.shape / 2
	for idx, v in ndenumerate(ft):
		ft[idx] = v * pow(1.j * (idx[0] - c[0]), dx) * pow(1.j * (idx[1] - c[1]), dy)
	return real(ifft2(ifftshift(ft)))

class Obstacle:

	def __init__(self, f, g, W, H, R):
		sf = array(imresize(f, W, H) > 0., 'f')
		sg = array(imresize(g, W, H) > 0., 'f')

		self.potential 	= se2_conv(f, g, R)
		grad_theta 		= se2_conv(diff_polar(f, 1, 0), g, R)

		grad_x 		= zeros((W, H, R))
		grad_y 		= zeros((W, H, R))
		for r in range(R):
			fr = imrotate(f, r / float(R) * 360.)
			fr_x = diff_cartesian(fr, 1, 0)
			fr_y = diff_cartesian(fr, 0, 1)
			grad_x[:,:,r] = real(fftconvolve(fr_x, g, 'same'))
			grad_y[:,:,r] = real(fftconvolve(fr_y, g, 'same'))

		self.grad = zeros((W, H, R, 3))
		self.grad[:,:,:,0] = grad_x
		self.grad[:,:,:,1] = grad_y
		self.grad[:,:,:,2] = grad_theta
		self.W = W
		self.H = H
		self.R = R

	def __xfrom_to_idx(self, config_f, config_g)
		rel = config_f.inv() * config_g
		if(abs(rel.x[0]) >= self.W/2 || abs(rel.x[1]) >= self.H/2 ):
			return None
		ix = rel.x.round()
		t = round(float(rel.theta) / (2. * pi) * float(self.R))
		t = (t + 10 * self.R) % self.R
		return ix, t

	def check_collide(self, config_f, config_g):
		return self.collision_potential(config_f, config_g) > 1.

	def collision_potential(self, config_f, config_g):
		ix, t = self.__xform_to_idx(config_f, config_g)
		if(ix is None):
			return 0
		return self.potential[ix[0], ix[1], t]

	def collision_gradient(self, config_f, config_g):
		ix, t = self.__xform_to_idx(config_f, config_g)
		if(ix is None):
			return 0
		return self.grad[ix[0], ix[1], t, :]

