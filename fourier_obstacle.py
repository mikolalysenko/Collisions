from scipy import array, exp, pi, sin, cos, arange, conjugate
import scipy.ndimage as ndi
from scipy.misc import imshow
import scipy.fftpack.pseudo_diffs as pds;
from precalc_obstacle import diff_polar
from polar import *
from se2 import *

#All shapes are fixed at the resolution, SHAPE_R, for the sake of simplicity
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
		
		self.mass_field = mass_field
		self.mass = sum(mass_field.flatten())

		#Recenter shape so that its center of mass is at the center of image
		center = array(ndi.center_of_mass(mass_field))
		offset = array(mass_field.shape)/2 - center
		nshape = (array(mass_field.shape) + 2. * abs(offset)).round()
		tmp = cpad(mass_field, nshape)
		self.mass_field = ndi.shift(tmp, offset, order=1)
		self.center = array(self.mass_field.shape) / 2
		
		assert(self.mass_field.shape[0] <= SHAPE_R)
		assert(self.mass_field.shape[1] <= SHAPE_R)
		
		self.indicator = array(self.mass_field > 0.01, 'f')
		
		#Compute moment of inertia
		self.moment = 0.
		self.radius = 0.
		for x,p in ndenumerate(self.mass_field):
			r = array(x) - self.center
			self.moment += p * dot(r,r)
			if(p > 0.01):
				self.radius = max(self.radius, norm(r))
				
		print self.radius
		self.shape_num = -1
		
		#Compute polar fourier truncation of indicator
		pind = cpad(self.indicator, array([2*SHAPE_R+1,2*SHAPE_R+1]))
		self.pft  = pfft(pind, R)
		self.pdft = map(pds.diff, self.pft)
		self.R = R
		
		#Compute residual energy terms
		self.energy = []
		s = real(self.pft[0][0]) * pi
		for r,l in enumerate(self.pft):
			s += sum(abs(real(self.pft[r]))) * r * 2. * pi / len(self.pft[r])
			self.energy.append(s)
		self.total_energy = s
		for r,e in enumerate(self.energy):
			self.energy[r] = self.total_energy - self.energy[r]

def c_shift(a, s):
	return pds.shift(a, s)

'''
The shape/obstacle data base
'''
class ShapeSet:
	def __init__(self, R):
		self.shape_list = []
		self.R = R

	'''
	Adds a shape to the obstacle set
	'''
	def add_shape(self, f):
		#Create shape
		S = Shape(f, self.R)
	
		#Add to shape list
		S.shape_num = len(self.shape_list)
		self.shape_list.append(S)

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
	

	'''
	Gradient calculation for shape field
	'''
	def grad(self, s1, s2, pa, pb, ra, rb):
		#Dereference shapes
		A = self.shape_list[s1]
		B = self.shape_list[s2]

		#Compute relative transformation
		ca = se2(pa, ra)
		cb = se2(pb, rb)
		rel = ca * cb.inv()
		
		if(max(abs(rel.x)) >= A.radius + B.radius ):
			return array([0.,0.,0.])
		
		#Load shape parameters
		fa  = A.pft
		fb  = B.pft
		db  = B.pdft
		ea 	= A.energy
		eb 	= B.energy
		cutoff = .01 * min(A.mass, B.mass) * (2. * SHAPE_R + 1) * (2. * SHAPE_R + 1)
		
		#Compute coordinate coefficients
		m   = 2.j * pi / (sqrt(2.) * SHAPE_R) * norm(rel.x)
		phi = atan2(rel.x[0], rel.x[1])
		
		#Set up initial sums
		s_0	= fa[0][0] * fb[0][0] * pi
		s_x = 0.
		s_y = 0.
		s_t = 0.
		
		for r in range(1, self.R):
		
			#Compute theta terms
			rscale = 2. * pi / len(fa[r])
			theta  = arange(len(fa[r])) * rscale
		
			#Construct multiplier / v
			mult =  fa[r] * exp((m * r) * cos(theta - phi)) * r * rscale
			v 	 =  pds.shift(fb[r], rel.theta) * mult
			
			#Check for early out
			s_0  += sum(real(v))
			if(s_0 + min(ea[r], eb[r]) <= cutoff):
				return array([0.,0.,0.])
				
			#Sum up gradient vectors
			v    *= 1.j
			s_x  += sum(real( v * cos(theta) ))
			s_y  -= sum(real( v * sin(theta) ))
			s_t	 -= r * sum(real(pds.shift(db[r], rel.theta) * mult))
		
		if(s_0 <= cutoff):
			return array([0., 0., 0.])
		return array([s_x, s_y, s_t])

