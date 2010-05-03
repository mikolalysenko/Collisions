'''
Misc. helper functions

-Mikola
'''
from scipy import array, zeros
from scipy.misc import imread

'''
Pads/centers an image
'''
def cpad(f, ns):
	res = zeros((ns[0],ns[1]), f.dtype)
	print res
	c0 = array((ns / 2. - array(f.shape)/2.).round(), 'i')
	c1 = c0 + array(f.shape, 'i')
	print c0, c1
	res[c0[0]:c1[0], c0[1]:c1[1]] = f
	return res
	
'''
Converts an image to an indicator map with the given cutoff
'''	
def to_ind(f, alpha=0.):
	return array(f > alpha, 'f')

'''
Loads an image from file
'''
def load_img(path):
	return to_ind(imread(path, flatten=True))


