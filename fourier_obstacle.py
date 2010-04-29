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
		self.pft  = pfft(cpad(self.indicator, array([2*SHAPE_R+1,2*SHAPE_R+1])), R)
		self.pdft = pfft(diff_polar(cpad(self.indicator, array([2*SHAPE_R+1,2*SHAPE_R+1])), 1), R)
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
		
	def potential(self, s1, s2, pa, pb, ra, rb):
		ca = se2(pa, ra)
		cb = se2(pb, rb)
		
		rel = ca * cb.inv()
		print rel
		A = self.shape_list[s1]
		B = self.shape_list[s2]
		if(ceil(abs(rel.x[0])) >= A.radius + B.radius or ceil(abs(rel.x[1])) >= A.radius + B.radius ):
			return 0.
		
		fa  = self.shape_list[s1].pft
		fb  = self.shape_list[s2].pft
		pr = norm(rel.x)
		phi = atan2(-rel.x[0], -rel.x[1])
		s_0	= fa[0][0] * fb[0][0] * pi
		m = 2.j * pi / (2*SHAPE_R + 1)
		for r in range(1, self.R):
			mult =  fa[r] * exp(m * pr * r * cos(arange(len(fa[r])) * 2. * pi / len(fa[r]) - phi))
			v 	 =  conjugate(shift(fb[r], rel.theta * len(fb[r]) / (2. * pi))) * mult
			s_0  += sum(v) * r * 2. * pi / (len(fb[r]))
		return real( s_0 ) / ((2. * SHAPE_R + 1.)**2)


	def grad(self, s1, s2, pa, pb, ra, rb):
		#if(s1 < s2):
		#	print "flip"
		#	return -self.grad(s2, s1, pb, pa, rb, ra)
		
		ca = se2(pa, ra)
		cb = se2(pb, rb)
		
		rel = ca * cb.inv()
		A = self.shape_list[s1]
		B = self.shape_list[s2]
		if(ceil(abs(rel.x[0])) >= A.radius + B.radius or ceil(abs(rel.x[1])) >= A.radius + B.radius ):
			return array([0.,0.,0.])
		
		print rel
		
		fa  = self.shape_list[s1].pft
		fb  = self.shape_list[s2].pft
		db  = self.shape_list[s2].pdft
		
		ea 	= self.shape_list[s1].energy
		eb 	= self.shape_list[s2].energy
		cutoff = .5 * self.shape_list[s1].total_energy * self.shape_list[s2].total_energy / (SHAPE_R * SHAPE_R)
		
		pr  = norm(rel.x)
		phi = atan2(rel.x[0], rel.x[1])
		s_0	= fa[0][0] * fb[0][0] * pi
		s_x = 0.
		s_y = 0.
		s_t = 0.
		m   = 1.j * pi / SHAPE_R
		for r in range(1, self.R):
			mult =  fa[r] * exp(m * pr * r * cos(arange(len(fa[r])) * 2. * pi / len(fa[r]) - phi))
			ss   =  0.
			v 	 =  conjugate(c_shift(fb[r], ss)) * mult
			h	 = 	r * 2. * pi / (len(fb[r]))
			s_0  += sum(real(v)) * h
			if(s_0 + ea[r] + eb[r] <= cutoff):
				print "early out",r,cutoff,s_0, ea[r], eb[r]
				return array([0.,0.,0.])
			s_x  += sum(real(v * 1.j * cos(arange(len(fb[r])) * 2. * pi / len(fb[r]))) ) * h
			s_y  += sum(real(v * 1.j * sin(arange(len(fb[r])) * 2. * pi / len(fb[r]))) ) * h
			s_t	 += sum(real(conjugate(c_shift(db[r], ss)) * mult))   * h
		
		#print s_0, s_x, s_y, s_t
		if(real(s_0) <= cutoff):
			return array([0., 0., 0.])
		g = array([real(s_x), real(s_y), real(s_t)])
		print "HIT: ", rel, g
		return -g
		
