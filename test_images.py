'''
Tests the obstacle calculation code

-Mikola
'''
from scipy import array, dtype, misc, pi, zeros, arange, cos, sin, real, imag, arange
from math import atan2
from scipy.misc import imshow, imsave
from rigid_body import *
from visualize import *
from misc import *

import fourier_obstacle as obstacle
from polar import ipfft, pft_mult
import os
import pickle

db = obstacle.ShapeSet(32, 256)
db.add_shape(load_img("shape3.png"))
db.add_shape(load_img("shape4.png"))


A = db.shape_list[0]
B = db.shape_list[1]

pconv = ipfft(pft_mult(A.pft, B.pft), 513, 513)
imsave("test1.png",  pconv)
imshow(pconv)

#Contstruct sampled convolution field
sconv = zeros((513, 513))
for i in range(513):
	for j in range(513):
		x = 257. - i
		y = 257. - j
		sconv[i, j] = db.potential(A, B, array([0., 0.]), array([x, y]), 0., 0.)
		print x,y, sconv[i,j]
		
imsave("test2.png", sconv)
imshow(sconv)

