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

db = obstacle.ShapeSet(32, 256)
db.add_shape(load_img("shape1.png"))
db.add_shape(load_img("shape2.png"))
db.add_shape(load_img("shape3.png"))
db.add_shape(load_img("shape4.png"))
db.add_shape(load_img("shape5.png"))

a = Body()
a.pos = array([-50.,130.])
a.rot = pi/2.
a.shape = db.shape_list[2]
a.lin_velocity = array([0., -40.])
a.ang_velocity = 0.0

b = Body()
b.pos = array([0.,0.])
b.rot = pi/2.
b.shape = db.shape_list[0]
b.lin_velocity = array([0., 0.])
b.ang_velocity = 0.

s = RigidBodySystem(db)
s.add_body(b)
s.add_body(a)

V = Visualization(s)
V.loop()

