import numpy as np
import matplotlib.pyplot as plt

ratio = [15,10,5,3,2,1.5,1.25,1.10]
top = [[34.282730942682974, 101.42900977110912, 111.97676378247888],
[34.999985868014846, 103.99058280668214, 114.61011797237916],
[37.27955890184967, 111.52220370724535, 120.816569977791],
[40.41695517823064, 120.3883968621958, 126.88621303728723],
[45.80214989681008, 131.80121184531345, 133.8221935710808],
[54.188484970944046, 145.87599711928488, 142.5387643609156],
[63.054748612110416, 159.21861113373598, 150.70205852937215],
[72.29431613732231, 172.3766754276403, 159.04948983549548],
]

bottom = [ [-192.98304260342957, -281.53077078634124, -299.38463459139894],
[-200.041743275302, -285.49241315053587, -306.10102720437317],
[-213.78373501571696, -295.37813619299413, -320.4098824817029],
[-227.3303754836379, -306.7341670262231, -335.78809513606575],
[-242.89838455375138, -319.9892549279157, -354.2115812569682],
[-259.5169802514465, -334.3066875346849, -374.5620308911093],
[-272.7892461422962, -344.81676661732126, -391.51649161963724],
[-284.55146882649206, -354.10410768667134, -407.1337225475145],
]


top_5 = [[-208.62999667939948, -148.2769066769597, -137.6497057990273],
[-273.1061676890131, -213.05097891015976, -201.31500678063276],
[-427.11094222036445, -366.5816808623222, -354.2063770351597],
[-596.1761655934606, -532.9462085927134, -520.4597710636262],
[-792.3074693760791, -725.5796998155329, -713.9296091187862],
[-992.7127283275586, -921.0413252549824, -910.6837740646479],
[-1162.0190639389257, -1085.8084527387946, -1078.1476423382155],
[-1314.1575885218735, -1232.9465640735034, -1227.4516038954498]
]

bottom_5 = [ [66.81802023195087, -28.156963199422535, -42.7310022844249],
[121.71655101828219, 29.844838810086003, 11.503009486496012],
[255.8321107184929, 170.366261541325, 144.62001511720518],
[401.81096137933673, 322.3332044917497, 290.5127008670752],
[571.8744979204876, 497.98611291469933, 459.9650462493719],
[744.4472900937653, 674.9888809876265, 631.7245866950143],
[892.2631141103684, 825.9794194758833, 777.7941318053818],
[-1314.1575885218735, -1232.9465640735034, -1227.4516038954498]]


# tmp = []
# for i,t in enumerate(top):
# 	tmp.append([abs(top_5[i][0]) - abs(t[0]),abs(top_5[i][1]) - abs(t[1]),abs(top_5[i][2]) - abs(t[2])]) 

tab_top = [ [abs(top_5[i][0]) - abs(t[0]),abs(top_5[i][1]) - abs(t[1]),abs(top_5[i][2]) - abs(t[2])] for i,t in enumerate(top)]

x = np.array([t[0] for t in tab_top])
x1 = np.array([t[1] for t in tab_top])
x2 = np.array([t[2] for t in tab_top])

# y = np.array(ratio)
# plt.plot(x, y,label="motor1")
# plt.plot(x1, y,label="motor2")
# plt.plot(x2, y,label="motor3")
# plt.xlabel("Forces (N)")
# plt.ylabel("Ratio Et/El")
# plt.legend()

# plt.show() # affiche la figure a l'ecran


youngModulusValues = [(94+i*216)*100/1390 for i in range(7)]
poissonRatioValues = [0.2616/(0.4161-i*0.02575) for i in range(7)]
angle_push_35_YM_PR = [56.8671,35.0064,21.154,13.1013,7.5752,3.2173,0.0053] # both YM & PR change
angle_pull_15_YM_PR = [22.2343,13.5283,9.1003,5.9182,3.4725,1.5672,0.000065669]
angle_pull_15_YM = [22.2343,13.8837,9.5698,6.5546,4.1945,2.1734,0.4247] # only YM change
angle_push_35_YM = [56.8671,36.1465,22.9118,14.9526,9.0325,4.194,0.0705]
angle_pull_15_PR = [22.2343,] # only PR change
angle_push_35_PR = [56.8671,56.2358,54.8738,55.1072,54.1273,53.0262,53.3378] 

# fig, axs = plt.subplots(1,2)
# fig.suptitle('Vertically stacked subplots')
# axs[0].plot(youngModulusValues, angle_push_35_YM_PR)
# axs[0].set_title('Axis [0, 0]')
# axs[1].plot(youngModulusValues, angle_pull_15_YM_PR)
# axs[1].set_title('Axis [0, 0]')

print(youngModulusValues)

from matplotlib.ticker import FormatStrFormatter
markers=['.',',','o','v','^','<','>','1','2','3','4','8','s','p','P','*','h','H','+','x','X','D','d','|','_']
descriptions=['point', 'pixel', 'circle', 'triangle_down', 'triangle_up','triangle_left', 'triangle_right', 'tri_down', 'tri_up', 'tri_left','tri_right', 'octagon', 'square', 'pentagon', 'plus (filled)','star', 'hexagon1', 'hexagon2', 'plus', 'x', 'x (filled)','diamond', 'thin_diamond', 'vline', 'hline']

fig, ax = plt.subplots()

# Set up grid, legend, and limits
ax.grid(True)
ax.plot(youngModulusValues, angle_push_35_YM_PR,label="pushing movement of 35 degree",linestyle='--', marker='o') #, color='b'
ax.plot(youngModulusValues, angle_pull_15_YM_PR,label="pulling movement of 15 degree",linestyle='--', marker='v') #\248
ax.set_xlabel("difference in pourcentage of El in conparison to Et")
ax.set_ylabel("total angle")
ax.legend()
ax.xaxis.set_major_formatter(FormatStrFormatter('%g'))
ax.xaxis.set_ticks(np.arange(0, 110, 10))

plt.show() # affiche la figure a l'ecran


# import matplotlib.pyplot as plt
# import numpy as np

# fig, ax = plt.subplots()

# ax.yaxis.set_major_formatter(FormatStrFormatter('%g'))
# ax.yaxis.set_ticks(np.arange(-2, 2, 0.25))

# x = np.arange(-1, 1, 0.1)
# plt.plot(x, x**2)
# plt.show()