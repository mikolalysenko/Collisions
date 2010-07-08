'''
Fourier obstacle for collision detection

-Mikola
'''
from scipy import array, exp, pi, sin, cos, arange, conjugate, real, sqrt, zeros
from math import atan2
from numpy import ndenumerate, dot
from scipy.linalg import norm
from scipy.misc import imshow, imrotate
from scipy.signal import fftconvolve

import scipy.ndimage as ndi
import scipy.fftpack.pseudo_diffs as pds


from misc import to_ind, cpad
from polar import pfft, ipfft, pft_mult, pft_rotate
from se2 import se2



'''
Computes the best cutoff for the given indicator function
'''
def best_cutoff(ift, pind, radius):
	pvalues = []
	for x,v in ndenumerate(ift):
		if(norm(array(x) - array(pind.shape) / 2.) <= radius):
			pvalues.append( (v, pind[x[0], x[1]]) )
	pvalues.sort()

	lmiss = zeros((len(pvalues)))
	umiss = zeros((len(pvalues)))
	
	l = 0
	u = 0
	for k in range(len(pvalues)):
		if(pvalues[len(pvalues) - k - 1][1] > 0):
			u += 1
		if(pvalues[k][1] == 0):
			l += 1
		lmiss[k] = l
		umiss[len(pvalues) - k - 1] = u
	descr = 3. * lmiss + umiss

	best_descr = descr[0]
	best_cutoff = 0.
	for k in range(len(pvalues)):
		if(descr[k] > best_descr):
			best_cutoff = pvalues[k][0]
			best_descr = descr[k]

	return best_cutoff


'''
Shape storage class
'''
class Shape:

	
	'''
	Initializes a shape object
	'''
	def __init__(self, mass_field, R, SHAPE_R):
		
		#Set basic parameters
		self.R = R
		self.SHAPE_R = SHAPE_R
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
		ift = real(ipfft(self.pft, pind.shape[0], pind.shape[1]))
		self.pft[0][0] -= min(ift.flatten()) * ((2. * SHAPE_R + 1) ** 2) #Enforce positivity
		self.pdft = map(pds.diff, self.pft)

		#Compute cutoff parameters
		ift = real(ipfft(self.pft, pind.shape[0], pind.shape[1]))
		self.cutoff = best_cutoff(ift, pind, self.radius)
		ind_ift = to_ind(ift, self.cutoff)
		self.int_res = sum((ift * ind_ift).flatten())
		self.ext_res = sum((ift * (1. - ind_ift)).flatten())
		self.res_area = sum(ind_ift.flatten())
		
		imshow(pind)
		imshow(to_ind(ift, self.cutoff))
		
		#Compute residual energy terms
		self.energy = []
		s = real(self.pft[0][0]) * pi
		for r,l in enumerate(self.pft):
			s += sum(abs(self.pft[r])) * (r * 2. * pi / len(self.pft[r]))
			self.energy.append(s)
		self.total_energy = s
		for r,e in enumerate(self.energy):
			self.energy[r] = s - self.energy[r]

'''
The shape/obstacle data base
'''
class ShapeSet:


	'''
	Initializes the shape database
	'''
	def __init__(self, R, SHAPE_R):
		self.shape_list = []
		self.cutoff_matrix = []
		self.R = R
		self.SHAPE_R = SHAPE_R
		self.tarea = ((2. * self.SHAPE_R + 1)**2)

	'''
	Adds a shape to the obstacle set
	'''
	def add_shape(self, f):
		#Create shape
		S = Shape(f, self.R, self.SHAPE_R)
	
		#Add to shape list
		S.shape_num = len(self.shape_list)
		self.shape_list.append(S)
		
		row = []
		for k in range(len(self.shape_list)):
			T = self.shape_list[k]
			ift = real(ipfft(pft_mult(pft_rotate(S.pft, 2.*pi/6.), T.pft), 2*self.SHAPE_R+1,2*self.SHAPE_R+1))
			Spad = imrotate(cpad(S.indicator, array([2*self.SHAPE_R+1,2*self.SHAPE_R+1])), 360./6.)
			Tpad = cpad(T.indicator, array([2*self.SHAPE_R+1,2*self.SHAPE_R+1]))
			pind = real(fftconvolve(Spad, Tpad, mode='same'))
			imshow(pind)
			imshow(ift)
			obst = to_ind(pind, 0.001)
			imshow(obst)
			cutoff = best_cutoff(ift, obst, S.radius + T.radius)
			print cutoff
			imshow(to_ind(ift, cutoff))
			row.append(cutoff * self.tarea)
		self.cutoff_matrix.append(row)

		return S

	
	'''
	Returns the cutoff threshold for shapes A,B
	'''
	def __get_cutoff(self, A, B):
		i = A.shape_num
		j = B.shape_num
		if(i < j):
			t = i
			i = j
			j = t
		C = self.cutoff_matrix[i][j]
		return C
		
	'''
	Evaluates shape potential field
	
	Not actually useful for collision detection, but somewhat helpful for debugging purposes.
	'''
	'''
	def potential(self, A, B, pa, pb, ra, rb):
		#Compute relative transformation
		ca = se2(pa, ra)
		cb = se2(pb, rb)
		rel = ca.inv() * cb
		pr = norm(rel.x)
		
		if(pr >= A.radius + B.radius ):
			return 0.
		
		#Load shape parameters
		fa  = A.pft
		fb  = B.pft
		ea 	= A.energy
		eb 	= B.energy
		cutoff = self.__get_cutoff(A,B)
		
		#Compute coordinate coefficients
		m   = 2.j * pi / (2. * self.SHAPE_R + 1) * pr
		phi = atan2(rel.x[1], rel.x[0])
				
		#Sum up energy contributions
		s_0	= real(fa[0][0] * fb[0][0] * pi)
		
		for r in range(1, self.R):
			#Compute theta terms
			rscale = 2. * pi / len(fa[r])
			theta  = arange(len(fa[r])) * rscale
		
			#Compute energy at this ring
			v =  pds.shift(fb[r], rel.theta) * fa[r] * exp((m * r) * cos(theta + phi)) * r * rscale
			
			#Check for early out
			s_0  += sum(real(v))
			if(s_0 + min(ea[r], eb[r]) <= cutoff):
				return 0.
		
		if(s_0 <= cutoff):
			return 0.
		return s_0
	'''
	
	'''
	Gradient calculation for shape field
	
	A,  B  - Shapes for the solids
	q      - Relative transformation
	'''
	def grad(self, A, B, q):
		#Compute relative transformation
		pr = norm(q.x)
		
		if(pr >= A.radius + B.radius ):
			return array([0., 0., 0., 0.])
		
		#Load shape parameters
		fa  = A.pft
		fb  = B.pft
		da  = A.pdft
		db  = B.pdft
		ea 	= A.energy
		eb 	= B.energy
		
		#Estimate cutoff threshold
		cutoff = self.__get_cutoff(A, B)
		
		#Compute coordinate coefficients
		m   = 2.j * pi / (2. * self.SHAPE_R + 1) * pr
		phi = atan2(q.x[1], q.x[0])
		
		#Set up initial sums
		s_0	 = real(fa[0][0] * fb[0][0])
		s_x  = 0.
		s_y  = 0.
		s_ta = 0.
		s_tb = 0.
		
		for r in range(1, self.R):
		
			#Compute theta terms
			dtheta = 2. * pi / len(fa[r])
			theta  = arange(len(fa[r])) * dtheta
		
			#Construct multiplier / v
			mult = exp((m * r) * cos(theta + phi)) * r * dtheta
			u 	 = pds.shift(conjugate(fb[r]), q.theta) * mult
			v 	 = fa[r] * u
			
			#Check for early out
			s_0  += sum(real(v))
			if(s_0 + min(ea[r], eb[r]) <= cutoff):
				return array([0.,0.,0.,0.])
				
			#Sum up gradient vectors
			v     = real(1.j * v)
			s_x  -= sum(v * sin(theta + phi) )
			s_y  -= sum(v * cos(theta + phi) )
			s_t  += sum(real(da[r] * u))
		
		if(s_0 <= cutoff):
			return array([0., 0., 0., 0.])
		return array([s_x, s_y, s_ta, s_tb, s_0])

