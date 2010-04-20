'''
Rigid body dynamics for python
'''
from scipy import matrix, array, eye
from obstacle import *
from se2 import *




class RigidBody:
	geometry = None;

	pos = None
	momentum = None
	
	def __init__(self, geom):
		assert(False);
	
	
