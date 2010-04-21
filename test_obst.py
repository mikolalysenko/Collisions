from scipy import array, dtype, misc, pi, zeros, arange, cos, sin, real, imag, arange
from scipy.misc import imshow

import obstacle

#Performs a level set (thresholding) operation on the given function
def to_ind(f, alpha=0.):
	return array(f > alpha, dtype('float'))

#Loads an image from file, converts it to an indicator function in R2
def load_img(path):
	return to_ind(misc.imread(path, flatten=True))

A = load_img("shape1.png")
S = shapes.add_shape(A)

print S.mass
print S.center
print S.moment

imshow(S.indicator)

import enthought.mayavi.mlab as mlab
mlab.pipeline.scalar_field(shapes.obstacles[0][0].potential)
mlab.pipeline.vector_field(shapes.obstacles[0][0].grad[:,:,:,0], shapes.obstacles[0][0].grad[:,:,:,1], shapes.obstacles[0][0].grad[:,:,:,2])



