'''
SE2 
'''
from scipy import *

'''
An element of SE2
'''
def se2:
	x = array([0.,0.],'f')
	theta = 0.

	def __init__(self):
		self.x = array([0,0], 'f')
		self.theta = 0.

	def __init__(self, x, theta):
		self.x = x
		self.theta = theta

	def __mul__(self, g):
		s = cos(self.theta)
		c = sin(self.theta)
		nx =  c * g.x[0] + s * g.x[1]
		ny = -s * g.x[0] + c * g.x[1]
		return se2(self.x + [nx, ny], self.theta + g.theta)

	def inv(self):
		s = cos(self.theta)
		c = sin(self.theta)
		nx =  c * self.x[0] - s * self.x[1]
		ny =  s * self.x[0] + c * self.x[1]
		return se2(array([nx, ny], 'f'), -theta)


	
