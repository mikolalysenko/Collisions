from scipy import array, dtype, misc, pi, zeros, arange, cos, sin, real, imag, arange
from math import atan2
from scipy.misc import imshow
from rigid_body import *
from visualize import *

import fourier_obstacle as obstacle

def to_ind(f, alpha=0.):
	return array(f > alpha, dtype('float'))

def load_img(path):
	return to_ind(misc.imread(path, flatten=True))

import os
import pickle

if(os.path.exists("db.pkl")):
	inp = open("db.pkl", "rb")
	db = pickle.load(inp)
	inp.close()
else:
	db = obstacle.ShapeSet(32)
	db.add_shape(load_img("shape1.png"))
	db.add_shape(load_img("shape2.png"))
	db.add_shape(load_img("shape3.png"))
	db.add_shape(load_img("shape4.png"))
	db.add_shape(load_img("shape5.png"))
	outp = open("db.pkl", "wb")
	pickle.dump(db, outp)
	outp.close()

print "here"

pp = zeros((64,64,3))
for x in range(64):
	for y in range(64):
		print x,y
		pp[x,y,:] = db.grad(4, 3, 8. * array([x,y],'f') - array([256., 256.]), array([0.,0.]), pi/2., 0.)
		print pp[x,y,:]
		#p = (array([x,y], 'f') - 32)
		#r = norm(p)
		#t = atan2(p[0], p[1])
		#pp[x,y] = eval_pft(s0.pft, r, t)
		#print pp[x,y]
print min(pp.flatten()), max(pp.flatten())
imshow(pp[:,:,2])

a = Body()
a.pos = array([10.,100.])
a.shape = db.get_shape(3)
a.lin_velocity = array([0., -40.])
a.ang_velocity = 0.0

b = Body()
b.pos = array([0.,0.])
b.rot = pi / 2.
b.shape = db.get_shape(4)
b.lin_velocity = array([0., 0.])
b.ang_velocity = 0.



s = RigidBodySystem(db)
s.add_body(a)
s.add_body(b)

print "here!"

V = Visualization(s)
V.loop()

