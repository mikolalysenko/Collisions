import shape;
from scipy import array, dtype, misc, pi;
from scipy.misc import imshow;

#Performs a level set (thresholding) operation on the given function
def to_ind(f, alpha=0.):
	return array(f > alpha, dtype('float'));

#Loads an image from file, converts it to an indicator function in R2
def load_img(path):
	return to_ind(misc.imread(path, flatten=True));


A = load_img("shape1.png");
imshow(A);

Apft = shape.geom2polar_ft(A, 10);
Atrunc = shape.polar_ft2geom(Apft, A.shape[0], A.shape[1]);
imshow(Atrunc);

for k in range(16):
	Arpft = shape.rotate_polar(Apft, pi * k / 32.);
	Artrunc = shape.polar_ft2geom(Arpft, A.shape[0], A.shape[1]);
	imshow(Artrunc);


