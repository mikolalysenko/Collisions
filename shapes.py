from scipy import *
from numpy import *
from scipy.ndimage import *
import obstacle

shape_list = []
obstacle_matrix = []

class Shape:

	def __init__(self, mass_field):
		self.mass_field = mass_field

		self.mass = sum(mass_field.flatten())
		center = center_of_mass(mass_field)

		#Recenter shape so that its center of mass is at the center of the coordinate system
		self.mass_field = shift(mass_field, (mass_field.shape - center) / 2)
		self.center = array(mass_field.shape) / 2

		#Compute moment of inertia
		self.moment = 0.
		for x,y,p in ndenumerate(self.mass_field):
			r = array([x,y]) - self.center
			moment += p * dot(r,r)

		self.indicator = array(mass_field > 0.001, 'f')

		#Add to shape list		
		self.shape_num = len(shape_list)
		shape_list.append(self)

		#Add new row to obstacle matrix
		obstacles = []
		for k,s in enumerate(shape_list):
			obstacles.append(Obstacle(self.indicator, s.indicator, R))
		obstacle_matrix.append(obstacles)

