# -*- coding: utf-8 -*-
import os
import math
from math import pi

from stlib.scene import Scene
from stlib.physics.constraints import FixedBox
from splib.numerics import *
import os, sys

from splib.numerics import sin, cos, to_radians
from stlib.physics.deformable import ElasticMaterialObject
from stlib.physics.collision import CollisionMesh

from splib.objectmodel import SofaPrefab, SofaObject
from splib.objectmodel import setData

from stlib.components import addOrientedBoxRoi
from splib.numerics import vec3
from splib.numerics.quat import Quat
from splib.numerics import *
from stlib.solver import DefaultSolver
from stlib.scene import Node

path = os.path.dirname(os.path.abspath(__file__))
path_utility = path+"/../utility/"
path_mesh = path + "/../mesh/"
sys.path.insert(0,path_utility)
from highOrderShape import HighOrderShape


# transverse_anisotropy_coef = [[4,1.0,0,0]]*nbr_tretra
# anisotropyDirections = [[0,1,0]]*nbr_tretra

# #######################################################################################################################################################
# # TRANSVERSE ONLY 

# # #FOR 15% AND SIGMA 0.1
# youngModulusTransverse = 1390
# poissonRatioTransverse = 0.26162 #(yz)
# youngModulusLongitudinal = 94.22
# poissonRatioTransverseLongitudinal = 0.4161  #(zx)
# shearModulusLongitudinal = 16.28215
# poissonRatioLongitudinalTransverse = 0.02975 #(xy)
# # ################################################

# # #FOR 15% AND SIGMA 0.1
# youngModulusTransverse = 1390
# poissonRatioTransverse = 0.26162 #(yz)
# youngModulusLongitudinal = 1390
# poissonRatioTransverseLongitudinal = 0.26162  #(zx)
# shearModulusLongitudinal = 550.879
# poissonRatioLongitudinalTransverse = 0.26162 #(xy)
# # ################################################


# youngModulus= [youngModulusTransverse]*nbr_tretra
# poissonRatio= [poissonRatioTransverse]*nbr_tretra

# for i,tetra in enumerate(directions_spiral):
#     anisotropyDirections[i] = directions_spiral[i]
#     transverse_anisotropy_coef[i] = [2,youngModulusLongitudinal,poissonRatioTransverseLongitudinal, shearModulusLongitudinal]
#     # transverse_anisotropy_coef[i] = [4,1,poissonRatioTransverseLongitudinal, shearModulusLongitudinal]

# #######################################################################################################################################################

# volumeMeshFileName = "mesh/cube_low_res.msh"
volumeMeshFileName = path_mesh + "beam_150-15-15.msh"
nbr_tretra = 189 #44


# print(transverse_anisotropy_coef)

path = os.path.dirname(os.path.abspath(__file__))
path_mesh = path + "/../mesh/"


youngModulus= None
poissonRatio= None
shearModulusLongitudinal = None
anisotropyDirections = None

transverse_anisotropy_coef = None


def createScene(rootNode):
    scene = Scene(rootNode, gravity=[0.0, 0.0,-9810],plugins=['SofaOpenglVisual','SoftRobots','SofaHighOrderTopology','SofaHighOrderFEM','BeamAdapter']) #,'SoftRobots.Inverse',"ConectToOptiTrack"])
    scene.VisualStyle.displayFlags = "showForceFields showBehavior"
    scene.createObject('BackgroundSetting', color=[0,0,0])

    totalMass=1

    # scene.createObject('FreeMotionAnimationLoop')
    # scene.createObject('GenericConstraintSolver')

    def cube(name="cube",translationY=0,rotation=[0,0,0],elasticitySymmetry='transverseIsotropic',ho=True,color=[0,0,0,1]):
        body = scene.createChild(name)

        if ho:

            e = HighOrderShape(body,name="body",
                            elasticitySymmetry=elasticitySymmetry,
                            translation=[0,-7.5+translationY,-7.5], rotation=rotation,
                            withAnysotropy=True,
                            volumeMeshFileName=volumeMeshFileName,
                            youngModulus=youngModulus, poissonRatio=poissonRatio, totalMass=totalMass,
                            anisotropyParameters=transverse_anisotropy_coef,anisotropyDirections = anisotropyDirections,color=[0,1,0],iterativeSolver=False)

            e.highOrder_node.createObject('EulerImplicitSolver', name='integration')
            e.highOrder_node.createObject('SparseLDLSolver', name="solver")
            # e.highOrder_node.createObject('GenericConstraintCorrection')
            e = e.highOrder_node



        else:
            e = ElasticMaterialObject(body,name="body",
                                    translation=[0,-7.5+translationY,-7.5], rotation=rotation,
                                    volumeMeshFileName=volumeMeshFileName,
                                    youngModulus=youngModulus, poissonRatio=poissonRatio, totalMass=totalMass)


        FixedBox(atPositions=[0, -8+translationY, -8, 1, 8+translationY, 8], applyTo=e,
                 doVisualization=True)

        setData(e.dofs, showObject=True, showObjectScale=1.8,drawMode=1, showColor=color)
 

    youngModulus= 1390
    poissonRatio= 0.26162    
    cube(name="isotropic",ho=False,color=[1,0,0,1])


    # transverseIsotropic # isotropic # cubic
    ####################################################
    youngModulus= [1390]*nbr_tretra
    poissonRatio= [0.26162]*nbr_tretra
    shearModulusLongitudinal = 550.879
    anisotropyDirections = [[1,1,0]]*nbr_tretra

    transverse_anisotropy_coef = [[4,1,poissonRatio[0], shearModulusLongitudinal]]*nbr_tretra
    cube(name="cubic",ho=True,elasticitySymmetry='cubic',color=[0,1,0,1])
    transverse_anisotropy_coef = [[2,youngModulus[0],poissonRatio[0], shearModulusLongitudinal]]*nbr_tretra
    cube(name="transverse",ho=True,elasticitySymmetry='transverseIsotropic',color=[0,0,1,1])
    ####################################################

    ####################################################
    youngModulus= 1390
    poissonRatio= 0.26162
    shearModulusLongitudinal = 550.879
    anisotropyDirections = [[1,0,0]]
    translationY=40

    transverse_anisotropy_coef = [4,1,poissonRatio, shearModulusLongitudinal]
    cube(name="iso",ho=True,elasticitySymmetry='isotropic',color=[1,0,0,1],translationY=translationY)
    cube(name="cubic",ho=True,elasticitySymmetry='cubic',color=[0,1,0,1],translationY=translationY)
    transverse_anisotropy_coef = [2,youngModulus-1e-11,poissonRatio, shearModulusLongitudinal] #-0.000000000001 # poissonRatio-0.0000000000000001
    cube(name="transverse",ho=True,elasticitySymmetry='transverseIsotropic',color=[0,0,1,1],translationY=translationY)
    ####################################################
