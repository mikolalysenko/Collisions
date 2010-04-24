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
S = obstacle.add_shape(A, 16)


print S.mass
print S.center
print S.moment

imshow(S.mass_field)
imshow(S.indicator)

O = obstacle.obstacle_matrix[0][0]
import enthought.mayavi.mlab as mlab
mlab.pipeline.scalar_field(O.potential)
mlab.pipeline.scalar_field(O.grad[:,:,:,2])
mlab.pipeline.scalar_field(O.grad[:,:,:,1])
mlab.pipeline.scalar_field(O.grad[:,:,:,0])
mlab.pipeline.vector_field(O.grad[:,:,:,0], O.grad[:,:,:,1], O.grad[:,:,:,2])


