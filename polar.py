from math import cos, sin, exp, sqrt, atan2, pi, floor;
from scipy import conjugate, real, array, dtype, zeros, concatenate, arange, ones;
from scipy.linalg.basic import norm;
from scipy.fftpack.basic import fft2, ifft2, fft, ifft;
from scipy.fftpack.helper import fftshift, ifftshift;
from scipy.fftpack.pseudo_diffs import shift;
from scipy.interpolate import interp1d;


#Non-uniform theta samples
__theta_samples = [ array([0]), array([0.*pi/8.,1.*pi/8.,2.*pi/8.,3.*pi/8.,4.*pi/8.,5.*pi/8.,6.*pi/8.,7.*pi/8.])]

'''
Interpolates non-uniform function f to circle
'''
def __resample_to_circle(f, R):
	s1 = concatenate((__theta_samples[R] - 2*pi, __theta_samples[R], __theta_samples[R] + 2*pi));
	ip = interp1d(s1,  concatenate((f, f, f)))
	res = zeros((8*R+4),f.dtype)
	for k in range(8*R + 4):
		theta = k * 2. * pi / float(8*R + 4);
		v = ip(theta);
		res[k] = v;
	return res;

'''
Resamples uniform f to a non-uniform function on the Bresenham circle
'''
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


'''
Initializes theta parameters for interpolation

FIXME: This should be precalculated somehow
'''
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
	x0 = int(round(.5 * xs))
	y0 = int(round(.5 * ys))
	
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
	
	x0 = int(round(.5 * xs))
	y0 = int(round(.5 * ys))

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



'''
Computes the (truncated) polar Fourier transform of f, with r rotational samples
'''
def pfft(f, nR):
	return rect2polar(ifftshift(fft2(ifftshift(__cpad(f, nR)))), nR);


'''
Computes the ivnerse polar Fourier transform of pft onto the xs by ys grid.
'''
def ipfft(pft, xs, ys):
	return fftshift(real(ifft2(fftshift(polar2rect(pft, xs, ys)))))

'''
Pads and centers f at a grid with sufficiently large rotational size
'''
def __cpad(f, r):
	ns = 2*max(int(round(norm(array(f.shape)/2.))), r)+1;
	res = zeros((ns,ns));
	c0 = array((ns / 2. - array(f.shape)/2.).round(), dtype('int'));
	c1 = c0 + array(f.shape, dtype('int'));
	res[c0[0]:c1[0], c0[1]:c1[1]] = f;
	return res;

'''
Performs pointwise polynomial multiplication to evaluate a single point
Point coordinates are given in polar form
'''
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
	
'''
Computes the dxth x derivative and the dyth y derivative of pft
'''
def pft_deriv(pft, dx, dy):
	res = []
	for r in range(len(pft)):
		c = pft[r].copy()
		for t in range(len(c)):
			theta = t * 2. * pi / len(c)
			c[t] *= pow(1.j * r * cos(theta), dx) * pow(1.j * r * sin(theta), dy)
		res.append(c)
	return res;
	
def pft_dpolar(pft, dr, dtheta):
	res = []
	for r in range(len(pft)):
		c = pft[r].copy()
		for t in range(len(c)):
			theta = t * 2. * pi / len(c)
			c[t] *= pow(1.j * r, dr) * pow(1.j * theta, dtheta)
		res.append(c)
	return res;
	
	
'''
Performs a polar convolution in the frequency domain
'''
def pft_mult(pX, pY):
	res = []
	for k in range(len(pX)):
		res.append(pX[k] * conjugate(pY[k]));
	return res;

'''
Rotates the uniform polar function pft
'''
def pft_rotate(pft, theta):
	R = len(pft);
	res = [];
	for k in range(R):
		res.append(shift(pft[k], theta));
	return res;

'''
Resamples f to resolution nr by interpolating in frequency
'''
def resample(f, nr):
	if(len(f) == 1):
		return ones(nr) * f[0];
	yy = fft(f)
	ft = interp1d(arange(len(yy)), yy, bounds_error=False, fill_value=0.)
	res = zeros((nr), yy.dtype)
	for k in range(nr):
		res[k] = ft(k)
	return ifft(res) * nr / len(f);



'''
Rescales the uniform polar function pft
'''
def pft_scale(pft, s):
	assert(s > 0);
	res = [];
	print len(pft);
	for k in range(int(len(pft) * s)):
		ns = int(round(k / s));
		if(ns >= len(pft)):
			break;
		res.append(__resample(pft[ns], 8*k+4));
	return res;

'''
Shifts the polar function pft
'''
def pft_shift(pft, x, y):
	assert(False);

'''
Performs scalar multiplication on pft
'''
def pft_scalar(pft, c):
	res = []
	for k in range(len(pft)):
		res.append(pft[k] * c);
	return res;

