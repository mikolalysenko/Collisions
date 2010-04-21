'''
Rigid body dynamics for python
'''
from scipy import matrix, array, eye
from shapes import *
from obstacle import *
from se2 import *



class RigidBody:
	shape_num = 0
	pos = None
	virtual_pos = None
	momentum = None
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
		return 0


