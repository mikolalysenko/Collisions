from scipy import array, dtype, misc, pi, zeros, arange, cos, sin, real, imag, arange
from scipy.misc import imshow

from obstacle import *
from rigid_body import *
from visualize import *

import pickle 
inp = open("shapes.pkl", "rb")
db = pickle.load(inp)
inp.close()

a = Body()
a.pos = array([400.,0.])
a.shape = db.get_shape(0)
a.lin_velocity = array([-1, 0])
a.ang_velocity = -0.001

b = Body()
b.pos = array([-400.,0.])
b.shape = db.get_shape(1)
b.lin_velocity = array([1, 0])
b.ang_velocity = 0.001

s = RigidBodySystem(db)
s.add_body(a)
s.add_body(b)

V = Visualization(s)
V.loop()

