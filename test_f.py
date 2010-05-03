'''
Test shell for setting up and running a basic system

-Mikola
'''
from scipy import array, dtype, misc, pi, zeros, arange, cos, sin, real, imag, arange
from math import atan2
from scipy.misc import imshow
from rigid_body import *
from visualize import *
from misc import *

import fourier_obstacle as obstacle
import polar as polar
import os
import pickle

if(os.path.exists("db.pkl")):
	inp = open("db.pkl", "rb")
	db = pickle.load(inp)
	inp.close()
else:
	db = obstacle.ShapeSet(32, 256)
	db.add_shape(load_img("shape1.png"))
	db.add_shape(load_img("shape2.png"))
	db.add_shape(load_img("shape3.png"))
	db.add_shape(load_img("shape4.png"))
	db.add_shape(load_img("shape5.png"))
	outp = open("db.pkl", "wb")
	pickle.dump(db, outp)
	outp.close()

a = Body()
a.pos = array([20.,100.])
a.shape = db.shape_list[3]
a.lin_velocity = array([0., -40.])
a.ang_velocity = 0.0

b = Body()
b.pos = array([0.,0.])
b.rot = pi / 2.
b.shape = db.shape_list[4]
b.lin_velocity = array([0., 0.])
b.ang_velocity = 0.

s = RigidBodySystem(db)
s.add_body(a)
s.add_body(b)

V = Visualization(s)
V.loop()

