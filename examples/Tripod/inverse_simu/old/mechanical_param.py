import numpy as np
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt

young_modulus_max = 1500
young_modulus_min = 90
young_modulus_step = 20
young_modulus_range = int((young_modulus_max - young_modulus_min) / young_modulus_step)

young_modulus_x = [] # [90+i for i in range(young_modulus_range)]
young_modulus_y = []
young_modulus_z = []

# print(young_modulus_range)
for i in range(young_modulus_range):
	for j in range(young_modulus_range):
		young_modulus_x.append(i*young_modulus_step+90)
		young_modulus_y.append(i*young_modulus_step+90)
		young_modulus_z.append(j*young_modulus_step+90)


# fig = plt.figure()
# ax = plt.axes(projection='3d')
# plt.title('Possible combinaison of young Modulus transverse')
# ax.set_xlabel('x')
# ax.set_ylabel('y')
# ax.set_zlabel('z')
# ax.plot3D(young_modulus_x, young_modulus_y, young_modulus_z, 'gray')
# plt.show()


poisson_ratio_max = 0.5
poisson_ratio_min = 0
poisson_ratio_step = 0.005
poisson_ratio_range = int((young_modulus_max - young_modulus_min) / young_modulus_step)

poisson_ratio_p = [] # [90+i for i in range(young_modulus_range)]
poisson_ratio_tp = []
poisson_ratio_pt = []

ratio = []
x = []
for i in range(young_modulus_range):
	x.append(i*young_modulus_step)
	ratio.append((young_modulus_min+i*young_modulus_step)/young_modulus_min)

plt.plot(x, ratio)
plt.show()

for i in range(poisson_ratio_range):
	for j in range(young_modulus_range):
		poisson_ratio_p.append(i*poisson_ratio_step+90)
		poisson_ratio_tp.append(i*young_modulus_step+90)
		poisson_ratio_pt.append(j*young_modulus_step+90)
