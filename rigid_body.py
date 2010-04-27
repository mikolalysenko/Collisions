'''
Rigid body dynamics for python
'''
from scipy import matrix, array, eye
from obstacle import *

class Body:
	pos = array([0., 0.])
	rot = 0.
	
	lin_velocity = array([0., 0.])
	ang_velocity = 0.
	
	force = array([0., 0.])
	torque = 0.

	shape = None

class RigidBodySystem:
	def __init__(self, shape_db):
		self.bodies = []
		self.gravity = array([0., -10.])
		self.shape_db = shape_db
	
	def add_body(self, body):
		self.bodies.append(body)
		
	def integrate(self, dt):
		for (i, body) in enumerate(self.bodies):
			body.v_pos = body.pos + dt * body.lin_velocity
			body.v_rot = body.rot + dt * body.ang_velocity
		
		for (i, A) in enumerate(self.bodies):
			for j in range(i):
				B = self.bodies[j]
				delta = self.shape_db.grad(A.shape.shape_num, B.shape.shape_num, A.v_pos, B.v_pos, A.v_rot, B.v_rot) * 100.
				if(abs(delta[0]) > 1):
					print i, j, delta
				A.force += delta[:2]
				A.torque += delta[2]
				B.force -= delta[:2]
				B.torque -= delta[2]
			A.force += self.gravity
			
		for (i, body) in enumerate(self.bodies):
			S = body.shape
			body.lin_velocity += body.force * dt / S.mass
			body.pos += body.lin_velocity * dt
			body.force = array([0.,0.])
			body.ang_velocity += body.torque * dt / S.moment
			body.rot += body.ang_velocity * dt
			body.torque = 0.
		
		return 0


