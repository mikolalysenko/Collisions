'''
Rigid body dynamics for python
'''
from scipy import matrix, array, eye
from shapes import *
from obstacle import *
from se2 import *



class RigidBody:
	shape_num = 0

	pos = array([0., 0.])
	rot = 0.
	
	lin_velocity = array([0., 0.])
	ang_velocity = 0.

	shape = None
	
	def __init__(self, geom):
		assert(False);


class RigidBodySystem:
	bodies = []
	gravity = -10.

	def __init__(self):
		self.bodies = []
		self.gravity = -10.
	
	def add_body(self, body):
		bodies.append(body)
	
	def integrate(self, dt):
		for (i, body) in enumerate(bodies):
			body.pos += dt * lin_velocity
			body.rot += dt * ang_velocity
		return 0


