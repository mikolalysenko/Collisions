'''
Rigid body dynamics for python
'''
from scipy import matrix, array, eye
from obstacle import *

class Body:
	shape_num = 0

	pos = array([0., 0.])
	rot = 0.
	
	lin_velocity = array([0., 0.])
	ang_velocity = 0.

	shape = None

class RigidBodySystem:
	bodies = []
	gravity = -10.

	def __init__(self, shape_db):
		self.bodies = []
		self.gravity = -10.
		self.shape_db = shape_db
	
	def add_body(self, body):
		self.bodies.append(body)
	
	def integrate(self, dt):
		for (i, body) in enumerate(self.bodies):
			body.pos += dt * body.lin_velocity
			body.rot += dt * body.ang_velocity
		return 0


