from os import *
from sys import *
from scipy import *
from numpy import *
import pygame
from pygame.locals import *
from OpenGL.GL import *
from OpenGL.GLU import *
import shapes
import obstacle
import rigid_body


class Visualization
	
	system = None
	textures = None

	def __init__(self, system):
		self.system = system
		self.screen = array([800,600])
		
		#Initialize pygame
		pygame.init()
		pygame.display.set_caption('Dynamics Visualizer')
		pygame.display.set_mode(self.screen,OPENGL|DOUBLEBUF)
		
		#Initialize opengl
		GL.resize(self.screen)
		GL.init()

		#Cache textures for all shapes

	def update(self):
		return 0

	def draw(self):
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
		pygame.display.flip()

	def loop(self):
		while(True):
			
			system.integrate(0.1)
			self.draw()
			


