shape_list = []
obstacle_matrix = []

class Shape:
	

	def __init__(self, mass_field):

		self.mass_field = mass_field
		self.center = array([0., 0.])
		self.moment = 1.
		self.indicator = array(mass_field > 0.001)

		
