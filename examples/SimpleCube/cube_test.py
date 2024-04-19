# -*- coding: utf-8 -*-
import os
import math
from math import pi

from stlib3.scene import Scene
from stlib3.physics.constraints import FixedBox
from splib3.numerics import *
import os, sys
from splib3.objectmodel import setData
from splib3.numerics import *
import Sofa

path = os.path.dirname(os.path.abspath(__file__))
path_mesh = path + "/mesh/"
# volumeMeshFileName = "mesh/cube_low_res.msh"
volumeMeshFileName = path_mesh + "beam_150-15-15.msh"


nbr_tretra = 189 #44

youngModulus= None
poissonRatio= None
shearModulusLongitudinal = None
anisotropyDirections = None

transverse_anisotropy_coef = None


def createScene(rootNode):

    def cube(name="cube",translationY=0,rotation=[0,0,0],elasticitySymmetry='transverseIsotropic',ho=True,color=[0,0,0,1]):
        # body = scene.createChild(name)

        if ho:
            e = Sofa.Core.Node(name)

            e.addObject('MeshGmshLoader',
                                      name='loader',
                                      rotation=rotation,
                                      translation=[0,-7.5+translationY,-7.5],
                                      filename=volumeMeshFileName)
            e.addObject('MeshTopology',
                                      src='@loader',
                                      name='container')

            e.addObject('MechanicalObject',
                                      name='dofs',
                                      position=e.loader.position.getLinkPath(),
                                      showObject=False,
                                      showObjectScale=5.0)
            e.addObject('UniformMass',
                                      name="mass",
                                      totalMass=0.032)

            # ForceField components

            high_fem = e.addObject('TetrahedronAnisotropicForceField',
                                      name="Elasticity",  printLog=1,  poissonRatio=poissonRatio,  youngModulus=youngModulus,
                                      # integrationMethod='analytical',
                                      # method='qr',
                                      # forceAffineAssemblyForAffineElements=False,
                                      # oneRotationPerIntegrationPoint=False,
                                      # numericalIntegrationMethod='Tetrahedron Gauss',
                                      # integrationOrder=1,
                                      # transparency= 0.75,
                                      drawHeterogeneousTetra=False)
            high_fem.elasticitySymmetry   = 'transverseIsotropic'
            high_fem.anisotropyParameters = transverse_anisotropy_coef
            high_fem.anisotropyDirections = anisotropyDirections



        else:
            e = Sofa.Core.Node(name)

            e.addObject('MeshGmshLoader',
                                      name='loader',
                                      rotation=rotation,
                                      translation=[0,-7.5+translationY,-7.5],
                                      filename=volumeMeshFileName)
            e.addObject('MeshTopology',
                                      src='@loader',
                                      name='container')

            e.addObject('MechanicalObject',
                                      name='dofs',
                                      position=e.loader.position.getLinkPath(),
                                      showObject=False,
                                      showObjectScale=5.0)
            e.addObject('UniformMass',
                                      name="mass",
                                      totalMass=0.032)

            # ForceField components

            e.addObject('TetrahedronFEMForceField',
                                      name="linearElasticBehavior",
                                      youngModulus=50,
                                      poissonRatio=0.45)
                                      # method="small")

        FixedBox(atPositions=[0, -8+translationY, -8, 1, 8+translationY, 8], applyTo=e,
                 doVisualization=True)

        setData(e.dofs, showObject=True, showObjectScale=1.8,drawMode=1, showColor=color)

        return e


    scene = Scene(rootNode, gravity=[0.0, 0.0,-9810], dt=0.001, iterative=True, plugins=["Anisotropy","SofaMatrix"])
    scene.addObject('BackgroundSetting', color=[1,1,1])
    scene.addMainHeader()
    scene.addObject('DefaultVisualManagerLoop')
    # scene.addObject('FreeMotionAnimationLoop')
    # scene.addObject('GenericConstraintSolver', maxIterations=50, tolerance=1e-5)
    # scene.Simulation.addObject('GenericConstraintCorrection')
    # scene.Simulation.TimeIntegrationSchema.rayleighStiffness = 0.005
    scene.Simulation.TimeIntegrationSchema.firstOrder=True
    # scene.Simulation.addObject("GlobalSystemMatrixImage")

    scene.Settings.mouseButton.stiffness = 10
    scene.VisualStyle.displayFlags = "showBehavior showForceFields"

    totalMass=1


    #### REF    ############################################
    youngModulus= 1390
    poissonRatio= 0.26162
    c = cube(name="isotropic",ho=False,color=[1,0,0,1])
    scene.Modelling.addChild(c)
    scene.Simulation.addChild(c)

    #### REF USING ANISOTROPY   ############################
    youngModulusTransverse = 1390 #1170
    poissonRatioTransverse = 0.26162 #(yz)
    youngModulusLongitudinal = 94.22
    poissonRatioTransverseLongitudinal = 0.4161  #(zx)
    shearModulusLongitudinal = 16.28215
    poissonRatioLongitudinalTransverse = 0.02975 #(xy)
    # youngModulusTransverse = 250
    # poissonRatioTransverse = 0.45 #(yz)
    # youngModulusLongitudinal = 250
    # poissonRatioTransverseLongitudinal = 0.45  #(zx)
    # shearModulusLongitudinal = 80
    # poissonRatioLongitudinalTransverse = 0.45 #(xy)
    ################################################

    youngModulus= [youngModulusTransverse]*nbr_tretra
    poissonRatio= [poissonRatioTransverse]*nbr_tretra
    anisotropyDirections = [[0.25,0.6,0]]*nbr_tretra
    transverse_anisotropy_coef = [[2,youngModulusLongitudinal,poissonRatioTransverseLongitudinal, shearModulusLongitudinal]]*nbr_tretra

    c = cube(name="isotropic - using plugin",ho=True,elasticitySymmetry='transverseIsotropic',color=[0,0,1,1])#,translationY=40)
    scene.Modelling.addChild(c)
    scene.Simulation.addChild(c)
    ####################################################


    # # transverseIsotropic # isotropic # cubic
    # ####################################################
    # youngModulus= [1390]*nbr_tretra
    # poissonRatio= [0.26162]*nbr_tretra
    # shearModulusLongitudinal = 550.879
    # anisotropyDirections = [[1,1,0]]*nbr_tretra
    #
    # transverse_anisotropy_coef = [[4,1,poissonRatio[0], shearModulusLongitudinal]]*nbr_tretra
    # cube(name="cubic",ho=True,elasticitySymmetry='cubic',color=[0,1,0,1])
    # transverse_anisotropy_coef = [[2,youngModulus[0],poissonRatio[0], shearModulusLongitudinal]]*nbr_tretra
    # cube(name="transverse",ho=True,elasticitySymmetry='transverseIsotropic',color=[0,0,1,1])
    # ####################################################
    #
    # ####################################################
    # youngModulus= 1390
    # poissonRatio= 0.26162
    # shearModulusLongitudinal = 550.879
    # anisotropyDirections = [[1,0,0]]
    # translationY=40
    #
    # transverse_anisotropy_coef = [4,1,poissonRatio, shearModulusLongitudinal]
    # cube(name="iso",ho=True,elasticitySymmetry='isotropic',color=[1,0,0,1],translationY=translationY)
    # cube(name="cubic",ho=True,elasticitySymmetry='cubic',color=[0,1,0,1],translationY=translationY)
    # transverse_anisotropy_coef = [2,youngModulus-1e-11,poissonRatio, shearModulusLongitudinal] #-0.000000000001 # poissonRatio-0.0000000000000001
    # cube(name="transverse",ho=True,elasticitySymmetry='transverseIsotropic',color=[0,0,1,1],translationY=translationY)
    # ####################################################
