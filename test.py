import polar;
import shape;
from scipy import array, dtype, misc, pi, zeros, arange, cos, sin, real, imag;
from scipy.misc import imshow;
from scipy.fftpack import fft2;
from scipy.signal.signaltools import convolve2d;


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
			res[xx,yy] = abs(shape.eval_pt(fpft, p * 0.1, phi));
			print res[xx,yy]
	return res;


A = load_img("shape1.png");
imshow(A);

B = load_img("shape2.png");
imshow(B);

Apft = shape.geom2polar_ft(A, 10);
Atrunc = shape.polar_ft2geom(Apft, A.shape[0], A.shape[1]);
imshow(Atrunc);

Bpft = shape.geom2polar_ft(B, 10);
Btrunc = shape.polar_ft2geom(Bpft, B.shape[0], B.shape[1]);
imshow(Btrunc);


'''
Asamp = plot_by_sample(Bpft, 64);
imshow(Asamp);
imshow(abs(Asamp));
imshow(real(Asamp));
imshow(imag(Asamp));
'''

C = convolve2d(A, B, mode='full');
imshow(C);

Cpft = shape.polar_conv(Apft, Bpft);
Ctrunc = shape.polar_ft2geom(Cpft, A.shape[0], A.shape[1]);
imshow(Ctrunc);

