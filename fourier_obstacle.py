from scipy import *
from numpy import *
import scipy.ndimage as ndi
from scipy.misc import *
from scipy.fftpack import *
from scipy.signal import *
from polar import *
from se2 import *


#All shapes are fixed at the resolution SHAPE_X, SHAPE_Y for simplicity
SHAPE_R = 256


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
	def __init__(self, mass_field, R):
		assert(mass_field.shape[0] == SHAPE_R)
		assert(mass_field.shape[1] == SHAPE_R)
			
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
		
		#Compute polar fourier truncation of signal
		self.pft = pfft_func(cpad(self.indicator, 2*SHAPE_R), R)


'''
The shape/obstacle data base
'''
class ShapeSet:
	def __init__(self):
		self.shape_list = []
		self.obstacle_matrix = []

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
		for T in self.shape_list:
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

