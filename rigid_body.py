'''
Rigid body dynamics for python
'''
from scipy import matrix, array, eye;
from polar import *;


'''
This class stores the geometric information for a specific rigid body
'''
class RigidGeom:
	mass = 1.;
	moment_of_inertia = 1.;
	pft_geom = [];
	
	'''
	Initializes rigid body geometry with given mass distribution
	'''
	def __init__(self, mass_field):
		assert(False);


class RigidBody:
	geometry = None;
	scale = 1.;
	
	position = array([0,0], dtype('float'));
	rotation = 0.;
	
	lin_momentum = array([0,0], dtype('float'));
	ang_momentum = 0.;
	
	def __init__(self, geom):
		assert(False);
	
	
	
