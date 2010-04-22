from scipy import array, dtype, misc, pi, zeros, arange, cos, sin, real, imag, arange
from scipy.misc import imshow

from obstacle import *
from rigid_body import *
from visualize import *

db = ShapeSet()
db.load_shapes("shapes.pkl")


b = Body()
b.shape = db.get_shape(0)
s = RigidBodySystem(db)
s.add_body(b)

V = Visualization(s)
V.loop()

