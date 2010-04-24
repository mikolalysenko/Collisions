from scipy import array, dtype, misc, pi, zeros, arange, cos, sin, real, imag, arange
from scipy.misc import imshow

import obstacle

#Performs a level set (thresholding) operation on the given function
def to_ind(f, alpha=0.):
	return array(f > alpha, dtype('float'))

#Loads an image from file, converts it to an indicator function in R2
def load_img(path):
	return to_ind(misc.imread(path, flatten=True))

db = obstacle.ShapeSet()

db.add_shape(load_img("shape1.png"), 64)
db.add_shape(load_img("shape2.png"), 64)

print len(db.shape_list)

import pickle 
outp = open("shapes.pkl", "wb")
pickle.dump(db, outp)
outp.close()
