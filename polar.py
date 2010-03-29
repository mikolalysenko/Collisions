from scipy import *

'''
Cartesian to Polar coordinate conversion

Uses Bresenham's algorithm to convert the image f into a polar sampled image
'''
def rect2polar( f, R ):
	#Check bounds on R
	assert(R > 0)

	xs, ys = f.shape
	x0 = floor(xs / 2)
	y0 = floor(ys / 2)
	
	fp = []
	
	#Initialize 0,1 as special cases
	fp.append(array([ f[x0, y0] ], f.dtype))
	
	if R >= 1:
		fp.append(array(
			[ f[x0,y0+1], f[x0+1,y0+1], f[x0+1,y0], f[x0+1,y0-1], f[x0,y0-1], f[x0-1,y0-1], f[x0-1,y0], f[x0-1,y0+1] ]))
			
	#Perform Bresenham interpolation
	for r in range(2, R):
		#Allocate result
		res = zeros((8 * r + 4), f.dtype)

		#Handle axial directions
		res[0+0*r] = f[x0,   y0+r]
		res[1+2*r] = f[x0+r, y0]
		res[2+4*r] = f[x0,   y0-r]
		res[3+6*r] = f[x0-r, y0]

		#Set up scan conversion process
		x = 0
		y = r
		s = 1 - r
		t = 1

		while x <= y:
			#Handle x-crossing
			x = x + 1
	
			res[    t+0] = f[x0+x,y0+y]
			res[2*r-t+1] = f[x0+y,y0+x]
			res[2*r+t+1] = f[x0+y,y0-x]
			res[4*r-t+2] = f[x0+x,y0-y]
			res[4*r+t+2] = f[x0-x,y0-y]
			res[6*r-t+3] = f[x0-y,y0-x]
			res[6*r+t+3] = f[x0-y,y0+x]
			res[8*r-t+4] = f[x0-x,y0+y]
			t = t + 1
	
			#Update status flag
			if  s < 0:
				s = s + 2 * x + 1
			elif x <= y:
				#Also handle y-crossing
				y = y - 1
				s = s + 2 * (x - y) + 1
	
				res[    t+0] = f[x0+x,y0+y]
				res[2*r-t+1] = f[x0+y,y0+x]
				res[2*r+t+1] = f[x0+y,y0-x]
				res[4*r-t+2] = f[x0+x,y0-y]
				res[4*r+t+2] = f[x0-x,y0-y]
				res[6*r-t+3] = f[x0-y,y0-x]
				res[6*r+t+3] = f[x0-y,y0+x]
				res[8*r-t+4] = f[x0-x,y0+y]
				t = t + 1
		fp.append(res)
        
	return fp


'''
Polar to Rectilinear coordinate conversion

Converts striped polar image into a Cartesian image
'''
def polar2rect( fp, xs, ys ):
	#Get size
	R = len(fp)
	assert(xs >= 2*R and ys >= 2*R)
	
	x0 = floor(xs/2)
	y0 = floor(ys/2)

	#Allocate result
	f = zeros((xs, ys), fp[0].dtype)
	
	#print(f.dtype)
	#print(fp[0].dtype)
	
	#Handle special cases
	f[x0, y0] = fp[0][0]

	if R >= 1:
		tmp = fp[1];
		f[x0,y0+1]   = tmp[0];
		f[x0+1,y0+1] = tmp[1];
		f[x0+1,y0]   = tmp[2];
		f[x0+1,y0-1] = tmp[3];
		f[x0,y0-1]   = tmp[4];
		f[x0-1,y0-1] = tmp[5];
		f[x0-1,y0]   = tmp[6];
		f[x0-1,y0+1] = tmp[7];

	#Perform scan conversion via Bresenham's algorithm
	for r in range(2, R):
		#Read circle values
		res = fp[r]
		
		#Set axial values
		f[x0,y0+r] = res[0+0*r]
		f[x0+r,y0] = res[1+2*r]
		f[x0,y0-r] = res[2+4*r]
		f[x0-r,y0] = res[3+6*r]
		
		#Begin Bresenham interpolation
		x = 0
		y = r
		s = 1 - r
		t = 1

		while x <= y:
			#Handle x-crossing
			x = x + 1
			
			f[x0+x,y0+y] = res[    t+0]
			f[x0+y,y0+x] = res[2*r-t+1]
			f[x0+y,y0-x] = res[2*r+t+1]
			f[x0+x,y0-y] = res[4*r-t+2]
			f[x0-x,y0-y] = res[4*r+t+2]
			f[x0-y,y0-x] = res[6*r-t+3]
			f[x0-y,y0+x] = res[6*r+t+3]
			f[x0-x,y0+y] = res[8*r-t+4]
			t = t + 1
			
			if s < 0:
				s = s + 2 * x + 1
			elif x <= y:
				#Also handle y-crossing
				y = y - 1
				s = s + 2 * (x - y) + 1
				f[x0+x,y0+y] = res[    t+0]
				f[x0+y,y0+x] = res[2*r-t+1]
				f[x0+y,y0-x] = res[2*r+t+1]
				f[x0+x,y0-y] = res[4*r-t+2]
				f[x0-x,y0-y] = res[4*r+t+2]
				f[x0-y,y0-x] = res[6*r-t+3]
				f[x0-y,y0+x] = res[6*r+t+3]
				f[x0-x,y0+y] = res[8*r-t+4]
				t = t + 1
	return f
