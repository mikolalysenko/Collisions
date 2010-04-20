'''
Rigid body dynamics for python
'''
from scipy import matrix, array, eye
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

	def add_body(self, body):
		bodies.append(body)
		
