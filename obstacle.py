import pickle
from scipy import *
from numpy import *
import scipy.ndimage as ndi
from scipy.misc import *
from scipy.fftpack import *
from scipy.signal import *
from polar import *
from se2 import *

'''
Implementation of SE2 convolution
'''
def se2_conv(f, g, R):
	r = np.zeros((f.shape[0], f.shape[1], R))
	for i in range(R):
		r[:,:,i] = real(fftconvolve(g, imrotate(f, i / float(R) * 360.), 'same'))
	return r

'''
Polar differentiation
'''
def diff_polar(fs, dtheta):
	R = int(round(sqrt(fs.shape[0] * fs.shape[0] + fs.shape[1] * fs.shape[1]) / 2 + 1.))
	f = cpad(fs, array([2*R+1,2*R+1]))
	pf = rect2polar(f, R)	
	for r in range(R):
		pf[r] = scipy.fftpack.diff(pf[r], dtheta)
	tmp = polar2rect(pf, 2*R+1, 2*R+1)
	tx = int(R - fs.shape[0]/2)
	ty = int(R - fs.shape[1]/2)
	return tmp[tx:(tx+fs.shape[0]), ty:(ty+fs.shape[1])]

'''
Cartesian differential
'''
def diff_cartesian(f, dx, dy):
	ft = fftshift(fft2(f))
	c = array(f.shape) / 2
	for idx, v in ndenumerate(ft):
		ft[idx] = v * pow(1.j * (idx[0] - c[0]), dx) * pow(1.j * (idx[1] - c[1]), dy)
	return real(ifft2(ifftshift(ft)))

'''
Pads/centers an image
'''
def cpad(f, ns):
	res = zeros((ns[0],ns[1]));
	c0 = array((ns / 2. - array(f.shape)/2.).round(), dtype('int'));
	c1 = c0 + array(f.shape, dtype('int'));
	res[c0[0]:c1[0], c0[1]:c1[1]] = f;
	return res;




'''
Shape storage class
'''
class Shape:
	def __init__(self, mass_field):
		self.mass_field = mass_field
		self.mass = sum(mass_field.flatten())

		#Recenter shape so that its center of mass is at the center of image
		center = array(ndi.center_of_mass(mass_field))
		offset = array(mass_field.shape)/2 - center
		nshape = (array(mass_field.shape) + 2. * abs(offset)).round()
		tmp = cpad(mass_field, nshape)
		self.mass_field = ndi.shift(tmp, offset, order=1)
		self.center = array(self.mass_field.shape) / 2

		#Compute moment of inertia
		self.moment = 0.
		for x,p in ndenumerate(self.mass_field):
			r = array(x) - self.center
			self.moment += p * dot(r,r)

		self.indicator = array(self.mass_field > 0.01, 'f')
		self.shape_num = -1

'''
Configuration obstacle / convolution field
'''
class Obstacle:
	def __init__(self, f, g, R):
		self.W = f.shape[0] + g.shape[1] + 1
		self.H = f.shape[1] + g.shape[1] + 1
		self.R = R

		sf = fliplr(flipud(cpad(f, array([self.W,self.H]))))
		sg = cpad(g, array([self.W,self.H]))

		self.potential 	= se2_conv(sf, sg, self.R)
		self.grad = zeros((self.W, self.H, self.R, 3))
		for r in range(self.R):
			print r
			fr = imrotate(sf, r / float(self.R) * 360.)
			fr_x = diff_cartesian(fr, 1, 0)
			fr_y = diff_cartesian(fr, 0, 1)
			fr_theta = diff_polar(fr, 1)
			self.grad[:,:,r,0] = real(fftconvolve(sg, fr_x, 'same'))
			self.grad[:,:,r,1] = real(fftconvolve(sg, fr_y, 'same'))
			self.grad[:,:,r,2] = real(fftconvolve(sg, fr_theta, 'same'))


	def xfrom_to_idx(self, config_f, config_g):
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
		ix, t = self.xform_to_idx(config_f, config_g)
		if(ix is None):
			return 0
		return self.potential[ix[0], ix[1], t]

	def collision_gradient(self, config_f, config_g):
		ix, t = self.xform_to_idx(config_f, config_g)
		if(ix is None):
			return 0
		return self.grad[ix[0], ix[1], t, :]

'''
The shape/obstacle data base
'''
class ShapeSet:
	shape_list = []
	obstacle_matrix = []

	'''
	Adds a shape to the obstacle set
	'''
	def add_shape(self, f, R):
		#Create shape
		S = Shape(f)
	
		#Add to shape list
		S.shape_num = len(self.shape_list)
		self.shape_list.append(S)

		#Generate obstacles
		obstacles = []
		for k,T in enumerate(self.shape_list):
			obstacles.append(Obstacle(S.indicator, T.indicator, R))
		self.obstacle_matrix.append(obstacles)
		return S
	
	'''
	Retrieves a shape with the given index
	'''
	def get_shape(self, idx):
		return self.shape_list[idx]
	
	def get_shapes(self):
		return self.shape_list

	def num_shapes(self):
		return len(self.shape_list)

	def grad(self, s1, s2, pa, pb, ra, rb):
		if(s1 < s2):
			return -self.check_collision(s2, s1, pb, pa, rb, ra)
		c1 = se2(pa, ra)
		c2 = se2(pb, rb)
		O = self.obstacle_matrix[s1][s2]
		return O.collision_gradient(c1, c2)

