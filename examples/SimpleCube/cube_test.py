# -*- coding: utf-8 -*-
import os, sys
import math

import Sofa
from stlib3.scene import Scene
from splib3.numerics import *
from splib3.numerics import *


from beam import *

import yaml
with open("mechaParam/anisotropyTransverse.yaml") as f:
    param = yaml.safe_load(f)

youngModulus = param["paramPaper"]["youngModulus"]["transverse"]

def createScene(rootNode):

    scene = Scene(rootNode, gravity=[0.0, 0.0,-98100], dt=0.01, iterative=True, 
      plugins=["Anisotropy","SofaMatrix"])
    # scene.addObject('BackgroundSetting', color=[1,1,1])
    scene.addMainHeader()
    scene.Simulation.TimeIntegrationSchema.rayleighMass = 0.1
    scene.Simulation.TimeIntegrationSchema.rayleighStiffness = 0.1
    scene.Simulation.TimeIntegrationSchema.vdamping = 1
    scene.addObject('DefaultVisualManagerLoop')
    # scene.addObject('FreeMotionAnimationLoop')
    # scene.addObject('GenericConstraintSolver', maxIterations=50, tolerance=1e-5)
    # scene.Simulation.addObject('GenericConstraintCorrection')
    # scene.Simulation.TimeIntegrationSchema.rayleighStiffness = 0.005
    # scene.Simulation.TimeIntegrationSchema.firstOrder=True
    # scene.Simulation.addObject("GlobalSystemMatrixImage")

    scene.Settings.mouseButton.stiffness = 10
    scene.VisualStyle.displayFlags = "showForceFields"



    toSimulate = []
    ################################################ REF Isotrope
    mechaParam = MechaParam(param["isotropic"])
    toSimulate.append(  Beam(name="isotropic",color=[1,0,0,1],
                        mechaParam=mechaParam,
                        elasticitySymmetry='isotropic').addModel()
                                                       .attach()  )

    ################################################ REF Isotrope with anisotropic ff
    mechaParam.setDirections([[0.25,0.5,0]])
    toSimulate.append(  Beam(name="isotropic",color=[1,0,0,1],
                        mechaParam=mechaParam,
                        elasticitySymmetry='transverseIsotropic').addModel()
                                                                 .attach()  )

    ################################################ REF Using Contol points
    mechaParam.setCP([[44,youngModulus, youngModulus ,radians(45),radians(0)],
                      [45,youngModulus, youngModulus ,radians(90),radians(0)]])
    toSimulate.append(  Beam(name="isotropic with CP",color=[1,0,0,1],
                        mechaParam=mechaParam,
                        elasticitySymmetry='transverseIsotropic').addModel()
                                                                 .attach()  )

    ################################################ Transverse right
    mechaParam = MechaParam(param["paramPaper"])
    mechaParam.setDirections([[0.25,0.5,0]])
    toSimulate.append(  Beam(name="Transverse right",color=[1,0,0,1],
                        mechaParam=mechaParam,
                        elasticitySymmetry='transverseIsotropic').addModel()
                                                                 .attach()  )

    ################################################ Transverse left
    mechaParam.setDirections([[-0.25,0.5,0]])
    toSimulate.append(  Beam(name="Transverse left",color=[1,0,0,1],
                        mechaParam=mechaParam,
                        elasticitySymmetry='transverseIsotropic').addModel()
                                                                 .attach()  )

    ################################################ Transverse Using Contol points
    mechaParam.setCP([[44,youngModulus, youngModulus-1200 ,radians(45),radians(0)],
                      [45,youngModulus-300, youngModulus-500,radians(90),radians(0)]])
    toSimulate.append(  Beam(name="Transverse CP 1",color=[1,0,0,1],
                        mechaParam=mechaParam,
                        elasticitySymmetry='transverseIsotropic').addModel()
                                                                 .attach()  )

    ################################################ Transverse Using Contol points
    mechaParam.setCP([[44,youngModulus, youngModulus-1200 ,radians(-45),radians(0)],
                      [45,youngModulus-300, youngModulus-500,radians(-90),radians(0)]])
    
    toSimulate.append(  Beam(name="Transverse CP 2",color=[1,0,0,1],
                        mechaParam=mechaParam,
                        elasticitySymmetry='transverseIsotropic').addModel()
                                                                 .attach()  )


    for beam in toSimulate: 
        scene.Modelling.addChild(beam)
        scene.Simulation.addChild(beam)
    