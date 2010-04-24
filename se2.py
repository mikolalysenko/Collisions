'''
SE2 
'''
from scipy import *

'''
An element of SE2
'''
class se2:
	def __init__(self):
		self.x = array([0,0], 'f')
		self.theta = 0.

	def __init__(self, x, theta):
		self.x = x
		self.theta = theta

	def __mul__(self, g):
		c = cos(self.theta)
		s = sin(self.theta)
		nx =  c * g.x[0] + s * g.x[1]
		ny = -s * g.x[0] + c * g.x[1]
		return se2(self.x + array([nx, ny]), self.theta + g.theta)

	def inv(self):
		c = cos(self.theta)
		s = sin(self.theta)
		nx =  c * self.x[0] - s * self.x[1]
		ny =  s * self.x[0] + c * self.x[1]
		return se2(array([-nx, -ny], 'f'), -self.theta)

	def __str__(self):
		return "[" + str(self.x) + "," + str(self.theta) + "]"

