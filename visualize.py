'''
A simple visualization template for rendering the results of the physics simulation

-Mikola
'''
from scipy import array, zeros, dtype, arange, pi
import pygame
from pygame.locals import *
from OpenGL.GL import *
from OpenGL.GLU import *

'''
The visualization system
'''
class Visualization:
	
	system = None
	textures = None
	running = True

	'''
	Initializes the visualization
	'''
	def __init__(self, system):
		self.system = system
		self.screen = array([800,800])
		
		#Initialize pygame
		pygame.init()
		pygame.display.set_caption('Dynamics Visualizer')
		pygame.display.set_mode(self.screen,OPENGL|DOUBLEBUF)
		
		glViewport(0, 0, self.screen[0], self.screen[1])

		glPixelStorei(GL_PACK_ALIGNMENT,1)			
		glPixelStorei(GL_UNPACK_ALIGNMENT,1)
		
		#HACK: PyOpenGL is stupid, sometimes it returns an array other times it doesn't.  WTF?
		shape_list = system.shape_db.shape_list
		num_shapes = len(shape_list)
		if(num_shapes == 0):
			self.textures = []
		elif(num_shapes == 1):
			self.textures = [glGenTextures(num_shapes)]
		else:
			self.textures = glGenTextures(num_shapes)
		
		#Cache all of the textures
		for s in shape_list:
			glBindTexture(GL_TEXTURE_RECTANGLE_ARB, self.textures[s.shape_num])
			glTexParameteri(GL_TEXTURE_RECTANGLE_ARB, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
			glTexParameteri(GL_TEXTURE_RECTANGLE_ARB, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
			glTexParameteri(GL_TEXTURE_RECTANGLE_ARB, GL_TEXTURE_WRAP_S, GL_CLAMP);
			glTexParameteri(GL_TEXTURE_RECTANGLE_ARB, GL_TEXTURE_WRAP_T, GL_CLAMP);
			
			#Need to flatten array and threshold colors properly
			s_flat = array(255. * s.indicator / max(s.indicator.flatten()), dtype('uint8'))
			glTexImage2D(GL_TEXTURE_RECTANGLE_ARB, 0, GL_LUMINANCE, s_flat.shape[1], s_flat.shape[0], 0, GL_RED, GL_UNSIGNED_BYTE, s_flat.tostring('C'))
			
	'''
	Handles user input events
	'''
	def do_events(self):
		for event in pygame.event.get():
			if event.type == QUIT:
				self.running = False
			elif event.type == KEYDOWN and event.key == K_ESCAPE:
				self.running = False
	
	'''
	Draws the simulation
	'''
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
		
		glDisable(GL_DEPTH_TEST)
		glEnable(GL_BLEND)
		glBlendFunc(GL_ONE, GL_ONE)
		
		glEnable(GL_TEXTURE_RECTANGLE_ARB)		
		for b in self.system.bodies:
			S = b.shape
			W = S.indicator.shape[1]
			H = S.indicator.shape[0]
			
			glBindTexture(GL_TEXTURE_RECTANGLE_ARB, self.textures[S.shape_num])
			glColor4f(1,1,1,1)
			
			glPushMatrix()
			glTranslatef(b.pos[0], b.pos[1], 0.)
			glRotatef(b.rot * 180. / pi, 0, 0, 1.)
			glScalef(.5, .5, 1.)
			
			glBegin(GL_QUADS)
			glTexCoord2f(0, 0)
			glVertex2f(-W, -H)
			
			glTexCoord2f(W, 0)
			glVertex2f( W, -H)
			
			glTexCoord2f(W, H)
			glVertex2f( W, H)
			
			glTexCoord2f(0, H)
			glVertex2f(-W, H)
			glEnd()
			
			glPopMatrix()
			
		glDisable(GL_BLEND)
		glDisable(GL_TEXTURE_RECTANGLE_ARB)
		pygame.display.flip()

	'''
	Main update loop
	'''
	def loop(self):
		while(self.running):
			self.do_events()
			self.system.integrate(0.001)
			self.draw()
			


