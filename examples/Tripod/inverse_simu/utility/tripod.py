import os, sys, math
from splib.numerics import sin, cos, to_radians
from stlib.physics.deformable import ElasticMaterialObject
from splib.objectmodel import setData

from actuatedarm import ActuatedArm

path = os.path.dirname(os.path.abspath(__file__))
path_utility = path+"/../../utility/"
path_mesh = path + "/../../mesh/"
sys.path.insert(0,path_utility)

from rigidification import Rigidify
from highOrderShape import HighOrderShape


def ElasticBody(parent,
				volumeMeshFileName= None,
				translation=[0.0,30,0.0], rotation=[90,0,0],
				youngModulus=600, poissonRatio=0.45, totalMass=0.4,**kwargs):

    body = parent.createChild("ElasticBody")

    e = ElasticMaterialObject(body,
                              volumeMeshFileName=volumeMeshFileName, # path_mesh + "tripod_coarse_04.vtk", #"tripod_mid.gidmsh",
                              translation=translation, rotation=rotation,
                              youngModulus=youngModulus, poissonRatio=poissonRatio, totalMass=totalMass)

    return body


def Tripod(parent, name="Tripod", radius=66, numMotors=3, angleShift=180.0, effectorPos=None, use_orientation=True, goalNode=None,ho=False,**kwargs):
    tripod = parent.createChild(name)

    maxForce = 150000 #20000
    # rest angle ~ 44degree choose range around rest of +/- 35degree
    maxPositiveDisp = 1.5708 #1.309 #1.5708  # 75degree
    maxNegativeDisp = 0 #-0.349066 #-0.226893 # 20degree #0
    maxDispVariation=0.2 #0.1
    initDisplacement = 0

    # maxPositiveDisp = 1.22173 #70
    # maxNegativeDisp = -0.349066 #20
    # maxDispVariation=0.2 #0.1
    # initDisplacement = 0

    # maxPositiveDisp = 1.0472 #60
    # maxNegativeDisp = -0.523599 #30
    # maxDispVariation= 0.2 #0.1
    # initDisplacement = 0

    # It is using the ElasticBody that was made in the previous step, and that
    # has also been included in the tripod.py script.
    if ho:
        body = tripod.createChild("ElasticBody")
        HighOrderShape(body,name="ElasticMaterialObject",**kwargs)
        body = tripod.ElasticBody.ElasticMaterialObject #.highOrder_node
    else:
        body = ElasticBody(tripod,**kwargs)
        body.init()
        body = tripod.ElasticBody.ElasticMaterialObject

    # The actuated arms are positionned around the silicone piece using a loop
    # structure
    dist = radius
    numstep = numMotors
    b = []
    arms = []
    angles = []
    frames = []
    test = []
    test1 = []
    for i in range(0, numstep):
        actuator_name = "ActuatedArm"+str(i)
        fi = float(i)
        fnumstep = float(numstep)
        angle = fi*360/fnumstep
        angle2 = fi*360/fnumstep+angleShift

        if name == "bottom":
	        ## BOTTOM
	        eulerRotation = [0, angle, 0]
	        angles.append([0, angle, 0])
	        translation = [dist*sin(to_radians(angle2)),
	                       -1.35-90,
	                       dist*cos(to_radians(angle2))]

	        frames.append([dist*sin(to_radians(angle2)),
	                       -1.35-90,
	                       dist*cos(to_radians(angle2)),
	                       0, angle, 0])
	        print(frames)
        elif name == "top":
	        ## TOP
	        eulerRotation = [0, -angle, 180]
	        angles.append([0, -angle, 0])
	        translation = [dist*sin(to_radians(angle2)),
	                       -1.35+90,
	                       dist*cos(to_radians(angle2))]

	        frames.append([dist*sin(to_radians(angle2)),
	                       -1.35+90,
	                       dist*cos(to_radians(angle2)),
	                       0, -angle, 180])
        else:
	        eulerRotation = [0, angle, ]
	        angles.append([0, angle, 0])
	        translation = [dist*sin(to_radians(angle2)),
	                       -1.35,
	                       dist*cos(to_radians(angle2))]

	        frames.append([dist*sin(to_radians(angle2)),
	                       -1.35,
	                       dist*cos(to_radians(angle2)),
	                       0, angle, ])

        #######################

        c = ActuatedArm(tripod, name=actuator_name,
                        translation=translation, eulerRotation=eulerRotation)

        arms.append(c)
        # if not ho:
        #     # print(body.dofs.rest_position)

        #     c.addBox(body.dofs.getData("rest_position"),
        #              translation, eulerRotation)
        #     # print(c.Box.BoxROI.indices)
        #     b.append(map(lambda x: x[0], c.Box.BoxROI.indices))
        # else:
        #     c.addBox(rest_position[i],
        #              translation, eulerRotation)
        #     b.append(map(lambda x: x[0], indices[i]))

    from math import pi
    if ho:
        if name == "bottom":
	        # frames = [[7.34788079488412e-15, -1.35-90, -60.0, 0, 0.0, 0], [-51.96152422706631, -1.35-90, 30.000000000000007, 0, 120.0, 0], [51.96152422706633, -1.35-90, 29.999999999999982, 0, 240.0, 0], [0, 30-90, 0, 0, 0, 0]]
	        frames = [[8.082668874372531e-15, -91.35, -66.0, 0, 0.0, 0], [-57.157676649772945, -91.35, 33.00000000000001, 0, 120.0, 0], [57.15767664977296, -91.35, 32.99999999999998, 0, 240.0, 0], [0, 30-90, 0, 0, 0, 0]] 
        	tmp = maxPositiveDisp
        	maxPositiveDisp = maxNegativeDisp
        	maxNegativeDisp = tmp
        	initDisplacement = -0.840425 # 48,15
        elif name == "top":
	        # frames = [[7.34788079488412e-15, -1.35+90, -60.0, 0, 0.0, 180], [-51.96152422706631, -1.35+90, 30.000000000000007, 0, -120.0, 180], [51.96152422706633, -1.35+90, 29.999999999999982, 0, -240.0, 180], [0, -30+90, 0, 0, 0, 180]]
	        frames = [[8.082668874372531e-15, 88.65, -66.0, 0, -0.0, 180], [-57.157676649772945, 88.65, 33.00000000000001, 0, -120.0, 180], [57.15767664977296, 88.65, 32.99999999999998, 0, -240.0, 180], [0, -30+90, 0, 0, 0, 180]]
        	initDisplacement = 0.784994 # 45
        else:
	        frames = [[7.34788079488412e-15, -1.35, -60.0, 0, 0.0, 0], [-51.96152422706631, -1.35, 30.000000000000007, 0, 120.0, 0], [51.96152422706633, -1.35, 29.999999999999982, 0, 240.0, 0], [0, 30, 0, 0, 0, 0]]
        
        b = [[14, 15, 16, 17, 100, 101, 102, 108, 109, 110, 116, 117, 118, 119, 120, 121, 122, 123, 129, 130, 131, 137, 138, 139, 222, 224, 225, 231, 232, 233, 234, 235, 239, 243, 244, 313, 349, 350, 362, 363, 388, 428, 429, 458, 461, 473, 476, 481, 483, 489, 493, 529, 565, 566, 578, 579, 604, 644, 645, 674, 677, 689, 692, 696, 697, 703, 707, 756, 760, 768, 784, 802, 818, 820, 825, 837, 842, 870, 874, 876, 878, 879], [2, 3, 4, 5, 25, 26, 27, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 49, 50, 51, 57, 58, 59, 176, 182, 183, 185, 186, 187, 188, 189, 190, 196, 197, 314, 348, 351, 360, 361, 389, 426, 427, 457, 459, 471, 475, 479, 480, 490, 491, 530, 564, 567, 576, 577, 605, 642, 643, 673, 675, 687, 691, 694, 695, 704, 705, 757, 759, 767, 785, 800, 807, 817, 821, 831, 836, 841, 843, 844, 859, 862, 871, 875, 877], [8, 9, 10, 11, 65, 66, 67, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 89, 90, 91, 97, 98, 99, 199, 206, 207, 208, 209, 210, 211, 212, 213, 219, 220, 312, 352, 353, 364, 365, 387, 425, 430, 456, 460, 472, 474, 478, 482, 488, 492, 528, 568, 569, 580, 581, 603, 641, 646, 659, 672, 676, 678, 688, 690, 702, 706, 758, 761, 783, 786, 801, 806, 808, 823, 824, 832, 857, 860, 880], [18, 19, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254, 255, 256, 257, 258, 259, 260, 261, 262, 263, 264, 265, 266, 267, 268, 269, 270, 271, 272, 273, 274, 275, 276, 277, 278, 279, 280, 281, 282, 283, 284, 285, 286, 287, 288, 289, 290, 291, 292, 293, 294, 295, 296, 297, 298, 299, 300, 301, 302, 378, 379, 380, 381, 382, 383, 384, 385, 386, 420, 421, 422, 423, 424, 431, 438, 439, 440, 737, 740, 741, 747, 749, 750, 751, 755, 773, 792, 794, 796, 797, 798, 799, 803, 804, 805, 809, 810, 811, 812, 813, 822, 827, 828, 840, 845, 846, 847, 848, 851, 852, 853, 856, 858, 861, 864, 866, 867, 868, 869]]


        o = Rigidify(tripod, body, groupIndices=b, frames=frames, name="RigidifiedStructure",noInit=True,ho=ho,pos=name)
    	setData(o.RigidParts.dofs, showObject=False, showObjectScale=10, drawMode=2)
    else:
        # if len(effectorPos) == 3:
        #     o = body.createObject("SphereROI", name="roi", template="Rigid3",
        #                 position=body.dofs.getData("rest_position"),
        #                 centers=effectorPos, radii=[7.5], drawSphere=True)
        #     o.init()
        #     b.append(map(lambda x: x[0], o.indices))

        #     frames.append([effectorPos[0],effectorPos[1],effectorPos[2],0,0,0])
            # print "----------------------------------"
            # print frames
            # print b
            # print "----------------------------------"
        if name == "bottom":
	        # frames = [[7.34788079488412e-15, -1.35-90, -60.0, 0, 0.0, 0], [-51.96152422706631, -1.35-90, 30.000000000000007, 0, 120.0, 0], [51.96152422706633, -1.35-90, 29.999999999999982, 0, 240.0, 0], [0, 30-90, 0, 0, 0, 0]]
	        frames = [[8.082668874372531e-15, -91.35, -66.0, 0, 0.0, 0], [-57.157676649772945, -91.35, 33.00000000000001, 0, 120.0, 0], [57.15767664977296, -91.35, 32.99999999999998, 0, 240.0, 0], [0, 30-90, 0, 0, 0, 0]] 
        elif name == "top":
	        # frames = [[7.34788079488412e-15, -1.35+90, -60.0, 0, 0.0, 180], [-51.96152422706631, -1.35+90, 30.000000000000007, 0, -120.0, 180], [51.96152422706633, -1.35+90, 29.999999999999982, 0, -240.0, 180], [0, -30+90, 0, 0, 0, 180]]
	        frames = [[8.082668874372531e-15, 88.65, -66.0, 0, -0.0, 180], [-57.157676649772945, 88.65, 33.00000000000001, 0, -120.0, 180], [57.15767664977296, 88.65, 32.99999999999998, 0, -240.0, 180], [0, -30+90, 0, 0, 0, 180]]
        else:
	        frames = [[7.34788079488412e-15, -1.35, -60.0, 0, 0.0, 0], [-51.96152422706631, -1.35, 30.000000000000007, 0, 120.0, 0], [51.96152422706633, -1.35, 29.999999999999982, 0, 240.0, 0], [0, 30, 0, 0, 0, 0]]
        
        b = [[14, 15, 16, 17, 100, 101, 102, 108, 109, 110, 116, 117, 118, 119, 120, 121, 122, 123, 129, 130, 131, 137, 138, 139, 222, 224, 225, 231, 232, 233, 234, 235, 239, 243, 244, 313, 349, 350, 362, 363, 388, 428, 429, 458, 461, 473, 476, 481, 483, 489, 493, 529, 565, 566, 578, 579, 604, 644, 645, 674, 677, 689, 692, 696, 697, 703, 707, 756, 760, 768, 784, 802, 818, 820, 825, 837, 842, 870, 874, 876, 878, 879], [2, 3, 4, 5, 25, 26, 27, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 49, 50, 51, 57, 58, 59, 176, 182, 183, 185, 186, 187, 188, 189, 190, 196, 197, 314, 348, 351, 360, 361, 389, 426, 427, 457, 459, 471, 475, 479, 480, 490, 491, 530, 564, 567, 576, 577, 605, 642, 643, 673, 675, 687, 691, 694, 695, 704, 705, 757, 759, 767, 785, 800, 807, 817, 821, 831, 836, 841, 843, 844, 859, 862, 871, 875, 877], [8, 9, 10, 11, 65, 66, 67, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 89, 90, 91, 97, 98, 99, 199, 206, 207, 208, 209, 210, 211, 212, 213, 219, 220, 312, 352, 353, 364, 365, 387, 425, 430, 456, 460, 472, 474, 478, 482, 488, 492, 528, 568, 569, 580, 581, 603, 641, 646, 659, 672, 676, 678, 688, 690, 702, 706, 758, 761, 783, 786, 801, 806, 808, 823, 824, 832, 857, 860, 880], [18, 19, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254, 255, 256, 257, 258, 259, 260, 261, 262, 263, 264, 265, 266, 267, 268, 269, 270, 271, 272, 273, 274, 275, 276, 277, 278, 279, 280, 281, 282, 283, 284, 285, 286, 287, 288, 289, 290, 291, 292, 293, 294, 295, 296, 297, 298, 299, 300, 301, 302, 378, 379, 380, 381, 382, 383, 384, 385, 386, 420, 421, 422, 423, 424, 431, 438, 439, 440, 737, 740, 741, 747, 749, 750, 751, 755, 773, 792, 794, 796, 797, 798, 799, 803, 804, 805, 809, 810, 811, 812, 813, 822, 827, 828, 840, 845, 846, 847, 848, 851, 852, 853, 856, 858, 861, 864, 866, 867, 868, 869]]


        o = Rigidify(tripod,
                    body,
                    name="RigidifiedStructure",
                    frames=frames,
                    groupIndices=b,noInit=True,ho=ho,pos=name)


# actuator0 = actuators.createObject('SlidingActuator', name='SlidingActuator0', template='Rigid3', direction='0 0 0 1 0 0', indices=0,
#                                         maxPositiveDisp=-0.5, maxNegativeDisp=1.5, maxDispVariation=0.02, initDisplacement=-1.47)
    actuators = o.RigidParts.createChild('actuators')
    actuator0 = actuators.createObject('SlidingActuator', name="SlidingActuator0", template='Rigid3d', 
    									direction='0 0 0 1 0 0' , indices=0, maxForce=maxForce, minForce=-maxForce, 
    									# maxDispVariation=maxDispVariation, initDisplacement=initDisplacement,
    									maxPositiveDisp=maxPositiveDisp,maxNegativeDisp=maxNegativeDisp)
    actuator1 = actuators.createObject('SlidingActuator', name="SlidingActuator1", template='Rigid3d',
    									direction='0 0 0 '+str(cos(4*math.pi/3))+' 0 '+str(sin(4*math.pi/3)) , indices=1, maxForce=maxForce, minForce=-maxForce, 
    									# maxDispVariation=maxDispVariation, initDisplacement=initDisplacement,
    									maxPositiveDisp=maxPositiveDisp,maxNegativeDisp=maxNegativeDisp)
    actuator2 = actuators.createObject('SlidingActuator', name="SlidingActuator2", template='Rigid3d',
    									direction='0 0 0 '+str(cos(2*math.pi/3))+' 0 '+str(sin(2*math.pi/3)) , indices=2, maxForce=maxForce, minForce=-maxForce,
    									# maxDispVariation=maxDispVariation, initDisplacement=initDisplacement,
    									maxPositiveDisp=maxPositiveDisp,maxNegativeDisp=maxNegativeDisp)

    # if goalNode==None:
    #     Effector = actuators.createObject('PositionEffector', name='effector', template='Rigid3d', useDirections='1 1 1 0 0 0', indices='3', effectorGoal="10 40 0", limitShiftToTarget=True, maxShiftToTarget=5 )
    # elif use_orientation:
    #     # Effector = actuators.createObject('PositionEffector', name='effector', template='Rigid3d', useDirections='0 1 0 1 0 1', indices='3', effectorGoal=goalNode.goalMO.getLinkPath()+".position", limitShiftToTarget=True, maxShiftToTarget=5)
    #     Effector = actuators.createObject('PositionEffector', name='effector', template='Rigid3d', useDirections='0 0 0 1 1 1', indices='3', effectorGoal=goalNode.goalMO.getLinkPath()+".position", limitShiftToTarget=True, maxShiftToTarget=5)
    # else:
    #     Effector = actuators.createObject('PositionEffector', name='effector', template='Rigid3d', useDirections='1 1 1 0 0 0', indices='3', effectorGoal=goalNode.goalMO.getLinkPath()+".position", limitShiftToTarget=True, maxShiftToTarget=5)

    # actuators.activated = 0


    sp = o.RigidParts.createChild("spring")
    sp.createObject('MechanicalObject', name='dofs', template='Vec3d',showObject=False,showObjectScale=20,
    	position=[[-18,0,0],[18,0,0],
				  [-18,0,0],[18,0,0],
				  [-18,0,0],[18,0,0]])
    sp.createObject("RigidMapping",rigidIndexPerPoint=[0,0,1,1,2,2])

    drawSpring = True
    sp.createObject('RestShapeSpringsForceField',
                              external_rest_shape=arms[0].ServoArm.spring.dofs.getLinkPath(),
							  points=[0,1],external_points=[0,1],
                              angularStiffness=1e6, stiffness=1e10,drawSpring=drawSpring)
    sp.createObject('RestShapeSpringsForceField',
                              external_rest_shape=arms[1].ServoArm.spring.dofs.getLinkPath(),
							  points=[2,3],external_points=[0,1],
                              angularStiffness=1e6, stiffness=1e10,drawSpring=drawSpring)
    sp.createObject('RestShapeSpringsForceField',
                              external_rest_shape=arms[2].ServoArm.spring.dofs.getLinkPath(),
							  points=[4,5],external_points=[0,1],
                              angularStiffness=1e6, stiffness=1e10,drawSpring=drawSpring)
    for i in range(0, numMotors):
        a = arms[i].ServoMotor.createChild("Attach")
        a.createObject("MechanicalObject", template="Rigid3d", name="dofs",
                       showObject=False, showObjectScale=10,
                       position=[0.0, 25.0, 10.0, 0, 0, 0, 1])
        a.createObject("RigidRigidMapping")

        # o.RigidParts.createObject('RestShapeSpringsForceField',
        #                           external_rest_shape=arms[i].ServoArm.dofs.getLinkPath(),
        #                           points=[i], external_points=[0],
        #                           angularStiffness=1e6, stiffness=1e10,drawSpring=True)

    for i in range(0, numMotors):
        arms[i].createObject("FixedConstraint")
        # arms[i].ServoMotor.ServoWheel.createObject("FixedConstraint")

    # tripod.removeChild(tripod.ElasticBody)

    return tripod