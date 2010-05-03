'''
Fourier obstacle for collision detection

'''
from scipy import array, exp, pi, sin, cos, arange, conjugate, real, sqrt
from math import atan2
from numpy import ndenumerate, dot
from scipy.linalg import norm

import scipy.ndimage as ndi
import scipy.fftpack.pseudo_diffs as pds

from misc import *
from polar import pfft
from se2 import se2

'''
Shape storage class
'''
class Shape:
	'''
	Initializes a shape object
	'''
	def __init__(self, mass_field, R, SHAPE_R):
		
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
		
		#Set indicator/shape area
		self.indicator = to_ind(self.mass_field, 0.01)
		self.area = sum(self.indicator.flatten())
		
		#Compute moment of inertia and radius
		self.moment = 0.
		self.radius = 0.
		for x,p in ndenumerate(self.mass_field):
			r = array(x) - self.center
			self.moment += p * dot(r,r)
			if(p > 0.01):
				self.radius = max(self.radius, norm(r))
		
		#Set shape indicator
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

'''
The shape/obstacle data base
'''
class ShapeSet:

	'''
	Initializes the shape database
	'''
	def __init__(self, R, SHAPE_R):
		self.shape_list = []
		self.R = R
		self.SHAPE_R = SHAPE_R

	'''
	Adds a shape to the obstacle set
	'''
	def add_shape(self, f):
		#Create shape
		S = Shape(f, self.R, self.SHAPE_R)
	
		#Add to shape list
		S.shape_num = len(self.shape_list)
		self.shape_list.append(S)

		return S

	
	'''
	Returns the cutoff threshold for shapes A,B
	'''
	def __get_cutoff(self, A, B):
		return .01 * min(A.area, B.area) * ((2. * self.SHAPE_R + 1) ** 2)
		
	'''
	Evaluates shape potential field
	'''
	'''
	def potential(self, A, B, pa, pb, ra, rb):
		#Compute relative transformation
		ca = se2(pa, ra)
		cb = se2(pb, rb)
		rel = ca * cb.inv()
		
		if(max(abs(rel.x)) >= A.radius + B.radius ):
			return 0.
		
		#Load shape parameters
		fa  = A.pft
		fb  = B.pft
		ea 	= A.energy
		eb 	= B.energy
		cutoff = self.__get_cutoff(A,B)
		
		#Compute coordinate coefficients
		m   = 2.j * pi / self.SHAPE_R * norm(rel.x)
		phi = atan2(rel.x[0], rel.x[1])
		
		#Sum up energy contributions
		s_0	= real(fa[0][0] * fb[0][0] * pi)
		
		for r in range(1, self.R):
			#Compute theta terms
			rscale = 2. * pi / len(fa[r])
			theta  = arange(len(fa[r])) * rscale
		
			#Compute energy at this ring
			v =  pds.shift(fb[r], rel.theta) * fa[r] * exp((m * r) * cos(theta - phi)) * r * rscale
			
			#Check for early out
			s_0  += sum(real(v))
			if(s_0 + min(ea[r], eb[r]) <= cutoff):
				return 0.
		
		if(s_0 <= cutoff):
			return 0.
		return 0.
	'''

	'''
	Gradient calculation for shape field
	
	A,  B  - Shapes for the solids
	pa, pb - Positions for shapes
	ra, rb - Rotations for shapes
	'''
	def grad(self, A, B, pa, pb, ra, rb):
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
		
		#Estimate cutoff threshold
		cutoff = self.__get_cutoff(A, B)
		
		#Compute coordinate coefficients
		m   = 2.j * pi / (sqrt(2.) self.SHAPE_R) * norm(rel.x)
		phi = atan2(rel.x[0], rel.x[1])
		
		#Set up initial sums
		s_0	= real(fa[0][0] * fb[0][0] * pi)
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

