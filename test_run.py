from scipy import array, dtype, misc, pi, zeros, arange, cos, sin, real, imag, arange
from scipy.misc import imshow
from rigid_body import *
from visualize import *

import pickle 
inp = open("shapes2.pkl", "rb")
db = pickle.load(inp)
inp.close()

for (i, s) in enumerate(db.shape_list):
	print i, s.shape_num


a = Body()
a.pos = array([0.,150.])
a.shape = db.get_shape(3)
a.lin_velocity = array([0., -10.])
a.ang_velocity = 0.0

b = Body()
b.pos = array([0.,0.])
b.rot = 0.
b.shape = db.get_shape(1)
b.lin_velocity = array([0., 0.])
b.ang_velocity = 0.

c = Body()
c.pos = array([0.,-250.])
c.shape = db.get_shape(0)
c.lin_velocity = array([0., 1.])
c.ang_velocity = 0.0

d = Body()
d.pos = array([200.,0])
d.shape = db.get_shape(1)
d.lin_velocity = array([-1., 0.])
d.ang_velocity = 0.0


s = RigidBodySystem(db)
s.add_body(a)
s.add_body(b)
#s.add_body(c)
#s.add_body(d)


V = Visualization(s)
V.loop()

