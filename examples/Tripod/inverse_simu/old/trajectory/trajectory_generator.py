
# Number_of_frame = 265
# Number_of_effector = 1
# KeyFrames = [0, 45, 92, 140, 187, 233, 264]

# ellipsoide = 

# trajectory = 
# "Type: Vec3d" + '\n' 
# "Number of frame: "    + str(Number_of_frame)	 + '\n'
# "Number of effector: " + str(Number_of_effector) + '\n'
# "KeyFrames: " 		   + str(KeyFrames)			 + '\n'

# for i in  


# Type: Vec3d
# Number of frame: 265
# Number of effector: 1
# KeyFrames: 0 45 92 140 187 233 264
# 0 0 170 
# -1.11111 0 170 
# -2.22222 0 170 
# -3.33333 0 170 
# -4.44444 0 170 

import numpy as np
# # define the lower and upper limits for x and y
# minX, maxX, minY, maxY, minZ, maxZ = 0., 120., 0., 80., 0., 80.
# # create one-dimensional arrays for x and y
# x = np.linspace(minX, maxX, 100) #(maxX-minX)/1000.+1)
# y = np.linspace(minY, maxY, 100) #(maxY-minY)/1000.+1)
# z = np.linspace(minZ, maxZ, 100) #(maxZ-minZ)/100.+1)

# # create the mesh based on these arrays
# X, Y, Z = np.meshgrid(x, y, z)
# coords = zip(X, Y, Z)
# print(X[0][0])
# print(Y[0][0])
# print(Z[0][0])

def coordWorkSpace()
	space = 100
	nbr_sample = 100/10 + 1
	x = [i*10 for i in range(nbr_sample)]
	coord = [[0,0,0]]*nbr_sample**3
	# print(space,nbr_sample)
	# print(x)
	count = 0
	for i in range(nbr_sample):
		for j in range(nbr_sample):
			for k in range(nbr_sample):
				# print(i,j,k)
				coord[count] = [x[i],x[j],x[k]]
				count +=1
	
	print(coord)
	return coord 