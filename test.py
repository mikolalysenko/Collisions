import polar;
from scipy import array, dtype, misc, pi, zeros, arange, cos, sin, real, imag, arange
from scipy.misc import imshow;
from scipy.fftpack import fft2;
from scipy.signal.signaltools import convolve2d;
from scipy.fftpack.pseudo_diffs import diff;

#Performs a level set (thresholding) operation on the given function
def to_ind(f, alpha=0.):
	return array(f > alpha, dtype('float'));

#Loads an image from file, converts it to an indicator function in R2
def load_img(path):
	return to_ind(misc.imread(path, flatten=True));

#Plots a function by samplin
def plot_by_sample(fpft, r):
	res = zeros((2*r+2, 2*r+2));
	for p in range(r):
		for phi in arange(0, 2*pi, 2*pi/(8*p+4)):
			xx = int(round(p * cos(phi) + r));
			yy = int(round(p * sin(phi) + r));
			print  "xx=", xx, "yy=",yy, "r=",p, "theta=",phi;
			res[xx,yy] = abs(polar.eval_pt(fpft, p * 0.1, phi));
			print res[xx,yy]
	return res;
	
	
import pylab;
x = arange(-pi, pi, 0.1);
y = sin(2. * x);


A = load_img("shape1.png");
imshow(A);

B = load_img("shape2.png");
imshow(B);

Apft = polar.pfft(A, 10);
Atrunc = polar.ipfft(Apft, A.shape[0], A.shape[1]);
imshow(Atrunc);

Bpft = polar.pfft(B, 10);
Btrunc = polar.ipfft(Bpft, B.shape[0], B.shape[1]);
imshow(Btrunc);

'''

Ashift = polar.pft_shift(Apft, 5., 0);
Astrunc = polar.ipfft(Ashift, A.shape[0], A.shape[1]);
imshow(Astrunc)

A2x = polar.pft_scale(Apft, 2.);
a2xt = polar.ipfft(A2x, A.shape[0], A.shape[1]);
imshow(a2xt);

Ahx = polar.pft_scale(Apft, .5);
ahxt = polar.ipfft(Ahx, A.shape[0], A.shape[1]);
imshow(ahxt);



for t in arange(0, pi/2., 0.1):
	imshow(polar.ipfft(polar.pft_rotate(Apft, t), A.shape[0], A.shape[1]));
'''

'''

Asamp = plot_by_sample(Bpft, 64);
imshow(Asamp);
imshow(abs(Asamp));
imshow(real(Asamp));
imshow(imag(Asamp));
'''


C = convolve2d(Atrunc, Btrunc);
imshow(C);

Cpft = polar.pft_mult(Apft, Bpft);
Ctrunc = polar.ipfft(Cpft, A.shape[0], A.shape[1]);
imshow(Ctrunc);

