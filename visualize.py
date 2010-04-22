from os import *
from sys import *
from scipy import *
from numpy import *
from ctypes import *
import pygame
from pygame.locals import *
from OpenGL.GL import *
from OpenGL.GLU import *
import obstacle
import rigid_body

class Visualization:
	
	system = None
	textures = None
	running = True

	def __init__(self, system):
		self.system = system
		self.screen = array([800,600])
		
		#Initialize pygame
		pygame.init()
		pygame.display.set_caption('Dynamics Visualizer')
		pygame.display.set_mode(self.screen,OPENGL|DOUBLEBUF)
		
		glViewport(0, 0, self.screen[0], self.screen[1])

		#Cache textures for all shapes
		glPixelStorei(GL_PACK_ALIGNMENT,1)			
		glPixelStorei(GL_UNPACK_ALIGNMENT,1)
		self.textures = glGenTextures(system.shape_db.num_shapes())
		for s in system.shape_db.get_shapes():
			glBindTexture(GL_TEXTURE_RECTANGLE_ARB, self.textures[s.shape_num])
			glTexParameteri(GL_TEXTURE_RECTANGLE_ARB, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
			glTexParameteri(GL_TEXTURE_RECTANGLE_ARB, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
			glTexParameteri(GL_TEXTURE_RECTANGLE_ARB, GL_TEXTURE_WRAP_S, GL_CLAMP);
			glTexParameteri(GL_TEXTURE_RECTANGLE_ARB, GL_TEXTURE_WRAP_T, GL_CLAMP);
			glTexImage2D(GL_TEXTURE_RECTANGLE_ARB, 0, GL_LUMINANCE, s.indicator.shape[0], s.indicator.shape[1], 0, GL_LUMINANCE, GL_UNSIGNED_BYTE, array(255. * s.indicator, dtype('byte')).flatten().tostring())

	def do_events(self):
		for event in pygame.event.get():
			if event.type == QUIT:
				self.running = False
			elif event.type == KEYDOWN and event.key == K_ESCAPE:
				self.running = False
	
	def draw(self):
		glClearColor(.23,.4,.8,0)
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
		
		glMatrixMode(GL_PROJECTION)
		glLoadIdentity()
		gluOrtho2D(-500, 500, -500, 500)
		
		glMatrixMode(GL_MODELVIEW)
		glLoadIdentity()
		
		glColor4f(1,1,1,1)
		glBegin(GL_LINES)
		for t in arange(-500., 500., 20.):
			glVertex2f(-500, t)
			glVertex2f( 500, t)
			glVertex2f(t, -500)
			glVertex2f(t,  500)
		glEnd()
		
		glEnable(GL_DEPTH_TEST)
		glEnable(GL_BLEND)
		glBlendFunc(GL_ONE, GL_ONE)
		
		glEnable(GL_TEXTURE_RECTANGLE_ARB)		
		for b in self.system.bodies:
			S = b.shape
			W = S.indicator.shape[0]
			H = S.indicator.shape[0]
			
			glBindTexture(GL_TEXTURE_RECTANGLE_ARB, self.textures[S.shape_num])
			glColor4f(1,1,1,1)
			
			glPushMatrix()
			glTranslatef(b.pos[0], b.pos[1], 0.)
			glRotatef(b.rot * 180. / pi, 0, 0, 1.)
			
			glBegin(GL_QUADS)
			glTexCoord2f(W, H)
			glVertex2f(-W/2., -H/2.)
			glTexCoord2f(0, H)
			glVertex2f( W/2., -H/2.)
			glTexCoord2f(W, 0)
			glVertex2f( W/2.,  H/2.)
			glTexCoord2f(0, 0)
			glVertex2f(-W/2.,  H/2.)
			glEnd()
			
			glPopMatrix()
			
		glDisable(GL_BLEND)
		glDisable(GL_TEXTURE_RECTANGLE_ARB)
		pygame.display.flip()

	def loop(self):
		while(self.running):
			self.do_events()
			self.system.integrate(0.1)
			self.draw()
			


