from math import cos, sin, exp, sqrt, atan2;
from scipy import conjugate, real, array;
from scipy.linalg.basic import norm;
from scipy.fftpack.basic import fft2, ifft2;
from scipy.fftpack.helper import fftshift, ifftshift;
from scipy.fftpack.pseudo_diffs import shift;
from scipy.interpolate import interp1d;

__theta_samples = [ array([0]), array([0.*pi/8.,1.*pi/8.,2.*pi/8.,3.*pi/8.,4.*pi/8.,5.*pi/8.,6.*pi/8.,7.*pi/8.])]

def __resample_to_circle(f, R):
	s1 = concatenate((__theta_samples[R] - 2*pi, __theta_samples[R], __theta_samples[R] + 2*pi));
	ip = interp1d(s1,  concatenate((f, f, f)))
	res = zeros((8*R+4),f.dtype)
	for k in range(8*R + 4):
		theta = k * 2. * pi / float(8*R + 4);
		v = ip(theta);
		res[k] = v;
	return res;

def __resample_from_circle(f, R):
	s = zeros((8*R+4))
	for k in range(8*R+4):
		s[k] = k * 2. * pi / float(8*R + 4);		
	s1 = concatenate((s-2*pi, s, s+2*pi))
	ip = interp1d(s1, concatenate((f, f, f)))
	res = zeros((8*R+4),f.dtype)
	for k in range(8*R+4):
		res[k] = ip(__theta_samples[R][k]);
	return res;


def __init_theta(nR):
	for r in range(len(__theta_samples),nR):
		#Allocate result
		res = zeros((8 * r + 4))
		p = 2 * r

		#Handle axial directions
		res[0    ] = 0.
		res[1+  p] = .5*pi
		res[2+2*p] = pi
		res[3+3*p] = 1.5*pi
		
		#Set up scan conversion process
		x = 0
		y = r
		s = 1 - r
		t = 1

		while t <= r:
			#Handle x-crossing
			x = x + 1
			
			theta = atan2(y, x)	
			res[    t+0] = res[1+  p] - theta
			res[  p-t+1] =              theta
			res[  p+t+1] = res[2+2*p] - theta
			res[2*p-t+2] = res[1+  p] + theta
			res[2*p+t+2] = res[3+3*p] - theta
			res[3*p-t+3] = res[2+2*p] + theta
			res[3*p+t+3] = 2*pi       - theta
			res[4*p-t+4] = res[3+3*p] + theta
			t = t + 1
	
			#Update status flag
			if  s < 0:
				s = s + 2 * x + 1
			elif t <= r:
				#Also handle y-crossing
				y = y - 1
				s = s + 2 * (x - y) + 1
				
				theta = atan2(y, x)	
				res[    t+0] = res[1+  p] - theta
				res[  p-t+1] =              theta
				res[  p+t+1] = res[2+2*p] - theta
				res[2*p-t+2] = res[1+  p] + theta
				res[2*p+t+2] = res[3+3*p] - theta
				res[3*p-t+3] = res[2+2*p] + theta
				res[3*p+t+3] = 2*pi       - theta
				res[4*p-t+4] = res[3+3*p] + theta
				t = t + 1
				
		print res;
		__theta_samples.append(res)
	


'''
Cartesian to Polar coordinate conversion

Uses Bresenham's algorithm to convert the image f into a polar sampled image
'''
def rect2polar( f, R ):
	#Check bounds on R
	assert(R > 0)
	
	__init_theta(R);

	xs, ys = f.shape
	x0 = floor(xs / 2)
	y0 = floor(ys / 2)
	
	fp = []
	
	#Initialize 0,1 as special cases
	fp.append(array([ f[x0, y0] ]))
	
	if R >= 1:
		fp.append(array(
			[ f[x0,y0+1], f[x0+1,y0+1], f[x0+1,y0], f[x0+1,y0-1], f[x0,y0-1], f[x0-1,y0-1], f[x0-1,y0], f[x0-1,y0+1] ]))
		
	#Perform Bresenham interpolation
	for r in range(2, R):
		#Allocate result
		res = zeros((8 * r + 4), f.dtype)
		p = 2 * r

		#Handle axial directions
		res[0    ] = f[x0,   y0+r]
		res[1+  p] = f[x0+r, y0]
		res[2+2*p] = f[x0,   y0-r]
		res[3+3*p] = f[x0-r, y0]

		#Set up scan conversion process
		x = 0
		y = r
		s = 1 - r
		t = 1

		while t <= r:
			#Handle x-crossing
			x = x + 1
	
			res[    t+0] = f[x0+x,y0+y]
			res[  p-t+1] = f[x0+y,y0+x]
			res[  p+t+1] = f[x0+y,y0-x]
			res[2*p-t+2] = f[x0+x,y0-y]
			res[2*p+t+2] = f[x0-x,y0-y]
			res[3*p-t+3] = f[x0-y,y0-x]
			res[3*p+t+3] = f[x0-y,y0+x]
			res[4*p-t+4] = f[x0-x,y0+y]
			t = t + 1
	
			#Update status flag
			if  s < 0:
				s = s + 2 * x + 1
			elif t <= r:
				#Also handle y-crossing
				y = y - 1
				s = s + 2 * (x - y) + 1
				
				res[    t+0] = f[x0+x,y0+y]
				res[  p-t+1] = f[x0+y,y0+x]
				res[  p+t+1] = f[x0+y,y0-x]
				res[2*p-t+2] = f[x0+x,y0-y]
				res[2*p+t+2] = f[x0-x,y0-y]
				res[3*p-t+3] = f[x0-y,y0-x]
				res[3*p+t+3] = f[x0-y,y0+x]
				res[4*p-t+4] = f[x0-x,y0+y]
				t = t + 1
				
		fp.append(__resample_to_circle(res, r))
        
	return fp


'''
Polar to Rectilinear coordinate conversion

Converts striped polar image into a Cartesian image
'''
def polar2rect( fp, xs, ys ):
	#Get size
	R = len(fp)
	assert(xs >= 2*R and ys >= 2*R)
	
	__init_theta(R);
	
	x0 = floor(xs/2)
	y0 = floor(ys/2)

	#Allocate result
	f = zeros((xs, ys), fp[0].dtype)
	
	#print(f.dtype)
	#print(fp[0].dtype)
	
	#Handle special cases
	f[x0, y0] = fp[0][0]


	if R >= 1:
		tmp = fp[1];
		f[x0,y0+1]   = tmp[0];
		f[x0+1,y0+1] = tmp[1];
		f[x0+1,y0]   = tmp[2];
		f[x0+1,y0-1] = tmp[3];
		f[x0,y0-1]   = tmp[4];
		f[x0-1,y0-1] = tmp[5];
		f[x0-1,y0]   = tmp[6];
		f[x0-1,y0+1] = tmp[7];
	
	#Perform scan conversion via Bresenham's algorithm
	for r in range(2, R):
		#Read circle values
		res = __resample_from_circle(fp[r], r)
		p = 2 * r
		
		#Set axial values
		f[x0,y0+r] = res[0    ]
		f[x0+r,y0] = res[1+  p]
		f[x0,y0-r] = res[2+2*p]
		f[x0-r,y0] = res[3+3*p]
		
		#Begin Bresenham interpolation
		x = 0
		y = r
		s = 1 - r
		t = 1

		while t <= r:
			#Handle x-crossing
			x = x + 1
			
			f[x0+x,y0+y] = res[    t+0]
			f[x0+y,y0+x] = res[  p-t+1]
			f[x0+y,y0-x] = res[  p+t+1]
			f[x0+x,y0-y] = res[2*p-t+2]
			f[x0-x,y0-y] = res[2*p+t+2]
			f[x0-y,y0-x] = res[3*p-t+3]
			f[x0-y,y0+x] = res[3*p+t+3]
			f[x0-x,y0+y] = res[4*p-t+4]
			t = t + 1
			
			if s < 0:
				s = s + 2 * x + 1
			elif t <= r:
				#Also handle y-crossing
				y = y - 1
				s = s + 2 * (x - y) + 1
				f[x0+x,y0+y] = res[    t+0]
				f[x0+y,y0+x] = res[  p-t+1]
				f[x0+y,y0-x] = res[  p+t+1]
				f[x0+x,y0-y] = res[2*p-t+2]
				f[x0-x,y0-y] = res[2*p+t+2]
				f[x0-y,y0-x] = res[3*p-t+3]
				f[x0-y,y0+x] = res[3*p+t+3]
				f[x0-x,y0+y] = res[4*p-t+4]
				t = t + 1
	return f



#Converts geometry into polar fourier transform
def pfft(real_geom, nR):
	return rect2polar(ifftshift(fft2(ifftshift(__cpad(real_geom)))), nR);

#Converts polar fourier transform into geometry
def ipfft(pft, xs, ys):
	return fftshift(real(ifft2(fftshift(polar2rect(pft, xs, ys)))))


#Pads geometry so that it doesn't have overlap issues
def __cpad(f):
	ns = int(round(norm(array(f.shape))));
	res = zeros((ns,ns));
	c0 = array((ns / 2. - array(f.shape)/2.).round(), dtype('int'));
	c1 = c0 + array(f.shape, dtype('int'));
	res[c0[0]:c1[0], c0[1]:c1[1]] = f;
	return res;

#Performs pointwise polynomial multiplication to evaluate a single point
# Point coordinates are given in polar form
def eval_pft(pft, p, phi):
	res = pft[0][0];	
	m = -2.j * pi / (len(pft));
	for r in range(1, len(pft)):
		c = pft[r];
		for t in range(len(c)):
			theta = t * 2*pi / len(c);
			v = exp(m * p * r * cos(theta - phi));
			#print cc, ss, v;
			res += c[t] * v;
	return res / (len(pft));
	
#Performs a convolution
def pft_mult(pX, pY):
	res = []
	for k in range(len(pX)):
		res.append(pX[k] * conjugate(pY[k]));
	return res;


#Rotates some polar geometry
def pft_rotate(pft, theta):
	R = len(pft);
	res = [];
	for k in range(R):
		res.append(shift(pft[k], theta));
	return res;

def pft_scale(pft, s):
	assert(False);
	
def pft_shift(pft, x, y):
	assert(False);


