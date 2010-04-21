from scipy import *
from numpy import *
from scipy.ndimage import *
from scipy.misc import *
from scipy.fftpack import *
from scipy.signal import *
from polar import *
from se2 import *

shape_list = []
obstacle_matrix = []



def se2_conv(f, g, R):
	r = np.zeros((f.shape[0], f.shape[1], R))
	for i in range(R):
		r[:,:,i] = real(fftconvolve(imrotate(f, i / float(R) * 360.), g, 'same'))
	return r

def diff_polar(f, dtheta):
	R = min(f.shape) / 2
	pf = rect2polar(f, R)	
	for r in range(R):
		pf[r] = scipy.fftpack.diff(pf[r], dtheta)
	tmp = polar2rect(pf, 2*R+1, 2*R+1)
	tx = int(R - f.shape[0]/2)
	ty = int(R - f.shape[1]/2)
	return tmp[tx:(tx+f.shape[0]), ty:(ty+f.shape[1])]

def diff_cartesian(f, dx, dy):
	ft = fftshift(fft2(f))
	c = f.shape / 2
	for idx, v in ndenumerate(ft):
		ft[idx] = v * pow(1.j * (idx[0] - c[0]), dx) * pow(1.j * (idx[1] - c[1]), dy)
	return real(ifft2(ifftshift(ft)))

def __cpad(f, ns):
	res = zeros((ns[0],ns[1]));
	c0 = array((ns / 2. - array(f.shape)/2.).round(), dtype('int'));
	c1 = c0 + array(f.shape, dtype('int'));
	res[c0[0]:c1[0], c0[1]:c1[1]] = f;
	return res;


class Shape:

	def __init__(self, mass_field):
		self.mass_field = mass_field
		self.mass = sum(mass_field.flatten())

		#Recenter shape so that its center of mass is at the center of the coordinate system
		offset = (array(mass_field.shape) - array(center_of_mass(mass_field))) / 2.
		self.mass_field = shift(mass_field, offset)
		self.center = array(mass_field.shape) / 2

		#Compute moment of inertia
		self.moment = 0.
		for x,p in ndenumerate(self.mass_field):
			r = array(x) - self.center
			self.moment += p * dot(r,r)

		self.indicator = array(mass_field > 0.001, 'f')
		self.shape_num = -1

class Obstacle:

	def __init__(self, f, g, R):
		self.W = f.shape[0] + g.shape[1] + 1
		self.H = f.shape[1] + g.shape[1] + 1
		self.R = R

		sf = fliplr(flipud(__cpad(f, array([W,H]))))
		sg = __cpad(g, array([W,H]))

		self.potential 	= se2_conv(f, g, R)
		self.grad = zeros((W, H, R, 3))
		for r in range(R):
			fr = imrotate(f, r / float(R) * 360.)
			fr_x = diff_cartesian(fr, 1, 0)
			fr_y = diff_cartesian(fr, 0, 1)
			fr_theta = diff_polar(fr, 1)
			self.grad[:,:,r,0] = real(fftconvolve(g, fr_x, 'same'))
			self.grad[:,:,r,1] = real(fftconvolve(g, fr_y, 'same'))
			self.grad[:,:,r,2] = real(fftconvolve(g, fr_theta, 'same'))


	def __xfrom_to_idx(self, config_f, config_g):
		rel = config_f.inv() * config_g
		if(abs(rel.x[0]) >= self.W/2 or abs(rel.x[1]) >= self.H/2 ):
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


'''
Adds a shape to the obstacle set
'''
def add_shape(f, R):
	S = Shape(f, R)

	S.shape_num = len(shape_list)
	shape_list.append(self)

	obstacles = []
	for k,T in enumerate(shape_list):
		obstacles.append(Obstacle(S.indicator, T.indicator, R))
	obstacle_matrix.append(obstacles)


