# -*- coding: utf-8 -*-
import os, sys
import math
from math import pi

import Sofa

from stlib3.scene import Scene
from splib3.numerics import *

from beam import *

import yaml
with open("mechaParam/anisotropyTransverse.yaml") as f:
    param = yaml.safe_load(f)

youngModulus = param["paramPaper"]["youngModulus"]["transverse"]


 

def createScene(rootNode):


    # rootNode.addObject('RequiredPlugin', pluginName=['SoftRobots','SoftRobots.Inverse','Anisotropy'])
    # rootNode.addObject('VisualStyle', displayFlags='showForceFields showBehaviorModels showCollisionModels')

    scene = Scene(rootNode, gravity=[0, -981.0, 0], dt=0.01, 
      plugins=['SoftRobots','SoftRobots.Inverse','Anisotropy'])
    scene.addMainHeader()
    scene.VisualStyle.displayFlags = "showForceFields showBehaviorModels showCollisionModels"

    rootNode.addObject('FreeMotionAnimationLoop',parallelODESolving=True)
    rootNode.addObject('QPInverseProblemSolver', printLog=False,epsilon=0)


    youngModulus= [1390]*nbr_tretra
    poissonRatio= [0.26162]*nbr_tretra
    shearModulusLongitudinal = 550.879
    anisotropyDirections = [[1,1,0]]*nbr_tretra

    transverse_anisotropy_coef = [[2,youngModulus[0],poissonRatio[0], shearModulusLongitudinal]]*nbr_tretra

    startingControlPoints=[[44,270, 80 ,radians(20),radians(60)],[51,500,  300 ,radians(10),radians(87)],[45,800,300,radians(31),radians(73)]]
    targetParamValue = [[44,700, 249 ,radians(90),radians(90)],[51,932,  500 ,radians(40),radians(80)],[45,462,212,radians(10),radians(90)]]

    ###
    ## WITH INVERSED BEAM
    ###
    masterActuator = "@/Simulation/c1_transverse/c1_transverse_model"
    anisotropyParameter = [[1,1,0,0]]*len(targetParamValue)
    toSimulate = []
    # for each cp need to give with the anisotropyParameter data a bool array indicating which value we will optimize:
    # the 2 first are for Y_t & Y_l the 2 other are for poissonRatioTransverse & poissonRatioTransverseLongitudinal
    # then we need to put the link (which will not be used) for all the slave AnisotropyControlPointActuator & put the one which will be "master"
    # at the beginning and it will go through the node and link with all the slaves during backwardInit()    
    ################################################ 
    posEffector=[[52.50000000000003, -7.5, 8.881784197001252e-16], [67.49999999999997, -7.5, -8.881784197001252e-16], [82.49999999999994, -7.5, 0.0], [97.49999999999989, -7.5, -2.6645352591003757e-15]]
    posGoal=[[48.490605069273535, -18.57535483905702, 0.6243967333346223], [62.152909669997, -24.10147114508542, 1.1016502954060348], [75.79754250445353, -29.883528746679747, 1.5811936429757942], [89.41867698686637, -35.984705877888366, 2.0197872173977633]]

    mechaParam = MechaParam(param["paramPaper"])
    mechaParam.setCP(startingControlPoints)
    toSimulate.append(  Beam(name="c1_transverse",color=[1,0,0,1],
                        mechaParam=mechaParam,translation=[0,-7.5,-7.5],
                        elasticitySymmetry='transverseIsotropic').addModel()
                                                                 .attach([0, -12, -12, 1, 12, 12])
                                                                 .addGoal(posGoal)
                                                                 .addEffector(posEffector)
                                                                 .addConstraint(targetParamValue,anisotropyParameter)
                                                                 )

    ################################################
    posEffector=[[52.50000000000003, -6.143643272632847, 20.698180332167443], [67.49999999999997, -6.143643272632844, 20.69818033216744], [82.49999999999994, -6.143643272632844, 20.69818033216744], [97.49999999999989, -6.143643272632843, 20.698180332167436]]
    posGoal=[[48.57398691126933, -18.10798009504002, 21.42377504374846], [62.09131226815256, -24.15813422510347, 21.769647797370038], [75.54888124335135, -30.470558588398998, 22.154120429565562], [88.9406452960323, -37.10410012527595, 22.5751696695389]]

    mechaParam = MechaParam(param["paramPaper"])
    mechaParam.setCP(startingControlPoints)
    toSimulate.append(  Beam(name="c2_transverse",color=[1,0,0,1],
                        mechaParam=mechaParam,translation=[0,-7.5+5.65818,-7.5-2.94546+25],rotation=[35,0,0],
                        elasticitySymmetry='transverseIsotropic').addModel()
                                                                 .attach([0, -12, -12+25, 1, 12, 12+25])
                                                                 .addGoal(posGoal)
                                                                 .addEffector(posEffector)
                                                                 .addConstraint(targetParamValue,anisotropyParameter,masterActuator=masterActuator) 
                                                                 )

    
    for beam in toSimulate: 
        scene.Modelling.addChild(beam)
        scene.Simulation.addChild(beam)
    

    ### WITH ONE RIGID AT TIP
    # posGoal=[[102.82297617948261, -102.19282738348612, 5.791339246658904, 0.026596082756856345, -0.04002919378724428, -0.4862751967142626, 0.8724830746650696]])

    #     c2 = cube(name="c2_transverse",ho=True,elasticitySymmetry='transverseIsotropic',color=[0,1,0,1],rotation=[35,0,0],tr=[0,5.65818,-2.94546],translationY=25,
    # ### WITH ONE RIGID AT TIP
    # # posGoal=[[97.4615348415232, -106.47374957061368, 31.993426182240015, 0.27007400681532295, -0.19650555226700042, -0.47600648936303896, 0.8135498762130737]])
