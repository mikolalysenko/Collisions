from scipy import array, dtype, misc, pi, zeros, arange, cos, sin, real, imag, arange
from scipy.misc import imshow
from rigid_body import *
from visualize import *

import pickle 
inp = open("shapes2.pkl", "rb")
db = pickle.load(inp)
inp.close()


a = Body()
a.pos = array([0.,150.])
a.shape = db.get_shape(3)
a.lin_velocity = array([0., -10.])
a.ang_velocity = 0.0

b = Body()
b.pos = array([0.,0.])
b.rot = pi / 2.
b.shape = db.get_shape(2)
b.lin_velocity = array([0., 0.])
b.ang_velocity = 0.

s = RigidBodySystem(db)
s.add_body(a)
s.add_body(b)

V = Visualization(s)
V.loop()

