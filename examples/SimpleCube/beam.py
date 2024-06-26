import os
import Sofa
from stlib3.physics.constraints import FixedBox
from splib3.objectmodel import setData
from splib3.numerics import *

rootPath = os.path.dirname(os.path.abspath(__file__))
path_mesh = rootPath + "/mesh/"

nbr_tretra = 189 #44

class MechaParam:
    def __init__(self,listParam,nbr_tretra=nbr_tretra):

        self.youngModulus = [listParam["youngModulus"]["transverse"]]*nbr_tretra
        self.youngModulusTransverse = listParam["youngModulus"]["transverse"]
        self.youngModulusLongitudinal = listParam["youngModulus"]["longitudinal"]

        self.poissonRatio = [listParam["poissonRatio"]["transverse"]]*nbr_tretra
        self.poissonRatioTransverse = listParam["poissonRatio"]["transverse"] #(yz)
        self.poissonRatioTransverseLongitudinal = listParam["poissonRatio"]["transverseLongitudinal"]  #(zx)
        self.poissonRatioLongitudinalTransverse = listParam["poissonRatio"]["longitudinalTransverse"] #(xy)
        
        self.shearModulusLongitudinal = listParam["shearModulus"]["longitudinal"]

        self.nbr_tretra = nbr_tretra

        self.transverse_anisotropy_coef = None

        self.setAnisotropyParam()

        self.anisotropyDirections = None
        self.controlPoints = None

    def setAnisotropyParam(self):
        self.transverse_anisotropy_coef = [[2,self.youngModulusLongitudinal,self.poissonRatioTransverseLongitudinal,self.shearModulusLongitudinal]]*self.nbr_tretra

    def setCP(self,controlPoints):
        self.controlPoints = controlPoints

    def setDirections(self,anisotropyDirections):
        self.anisotropyDirections = anisotropyDirections*self.nbr_tretra

class Controller(Sofa.Core.Controller):

    def __init__(self, *args, **kwargs):
        # These are needed (and the normal way to override from a python class)
        Sofa.Core.Controller.__init__(self, *args, **kwargs)

        self.posIndices = [52,53,49,47]
        
        self.goalMecha = kwargs["goal"]
        self.posGoal = None
        self.posEffector = None

    def onKeypressedEvent(self, event):
        key = event['key']
        if ord(key) == 19:  # up
            pos = self.getContext().dofs.position.value
            self.posEffector=[list(pos[i]) for i in self.posIndices]
            self.posGoal = [list(i) for i in self.goalMecha.position.value]

            print("###########################################")
            print(self.name.value)
            print("posEffector="+str(self.posEffector))
            print("posGoal="+str(self.posGoal))
            print("delta="+str([[l-k for l,k in zip(i,j)] for i,j in zip(self.posEffector,self.posGoal)]))
            print("###########################################")


def Beam(mechaParam,parent=None,name='Beam',mesh="beam_150-15-15.msh",elasticitySymmetry='transverseIsotropic',
         rotation=[0,0,0],scale3d=[1,1,1],translation=[0,0,0],color=[0,0,0,1]):

    if (parent == None):
        self = Sofa.Core.Node(name)
    else:
        self = parent.addChild(name)

    self.rotation = rotation
    self.translation = translation
    self.scale3d = scale3d

    self.elasticitySymmetry = elasticitySymmetry
    self.mechaParam = mechaParam

    self.addObject('VisualStyle',name='visualStyle')
    
    def addModel():


        self.model = self.addChild(self.name.value+"_model")


        self.model.addObject('MeshGmshLoader',
                                  name='loader',
                                  rotation=rotation,
                                  translation=translation,
                                  filename=path_mesh+mesh)
        self.model.addObject('MeshTopology',
                                  src='@loader',
                                  name='container')

        self.model.addObject('MechanicalObject',
                                  name='dofs',
                                  position=self.model.loader.position.getLinkPath())

        self.model.addObject('UniformMass',
                                  name="mass",
                                  totalMass=1)#0.032)
        if self.elasticitySymmetry == "transverseIsotropic":


            # ForceField components
            if (mechaParam.controlPoints):
                ff = self.model.addObject('TetrahedronAnisotropicForceField', printLog=1,
                                        name="Elasticity", 
                                        elasticitySymmetry= 'transverseIsotropic',
                                        controlPoints=mechaParam.controlPoints,
                                        meshRotation=radians(self.rotation[0]),
                                        drawHeterogeneousTetra=True)
            else: 
                ff = self.model.addObject('TetrahedronAnisotropicForceField', printLog=1,
                                        name="Elasticity",  
                                        poissonRatio=mechaParam.poissonRatio, 
                                        youngModulus=mechaParam.youngModulus,
                                        drawHeterogeneousTetra=False)
                print(ff,type(ff))
                ff.elasticitySymmetry   = 'transverseIsotropic'
                ff.anisotropyParameters = mechaParam.transverse_anisotropy_coef
                ff.anisotropyDirections = mechaParam.anisotropyDirections

        else:

            self.model.addObject('TetrahedronFEMForceField',#printLog=1
                                      name="linearElasticBehavior",
                                      youngModulus=mechaParam.youngModulus,
                                      poissonRatio=mechaParam.poissonRatio,
                                      method="large")

        # setData(self.model.dofs, showObject=True, showObjectScale=1,drawMode=1, showColor=color)
        
        return self

    def attach(atPositions=[0, 0, 0, 1, 16, 16]):

        FixedBox(atPositions=atPositions, applyTo=self.model,
                 doVisualization=True)
        return self

    def addGoal(posGoal):

        self.goal = self.addChild(self.name.value+'_goal')
        self.goal.addObject('EulerImplicitSolver', firstOrder=True)
        self.goal.addObject('CGLinearSolver', iterations='100', tolerance="1e-5", threshold="1e-5")

        ### WITH VEC3
        self.goal.addObject('MechanicalObject', name='goalMO',position=posGoal)#[posGoal[0][:3]])
        self.goal.addObject('SphereCollisionModel', radius=1, group='3')

        ### WITH ONE RIGID AT TIP
        # goal.addObject('MechanicalObject', name='goalMO',template='Rigid3d',
        #                   position=posGoal, showObject=True, showObjectScale=3,drawMode=3)
        ###
        self.goal.addObject('UncoupledConstraintCorrection')   
        return self

    def addEffector(posEffector):

        effector = self.model.addChild('effector')
        effector.addObject('MechanicalObject', name="effectorPoint",showObject=True,showObjectScale=10,position=posEffector)
        effector.addObject('PositionEffector', template='Vec3', indices=[0,1,2,3], effectorGoal=self.goal.goalMO.getLinkPath()+".position")
        effector.addObject('BarycentricMapping', mapForces=False, mapMasses=False)

        return self

    def addConstraint(targetParamValue,anisotropyParameter=[[1,1,0,0],[1,1,0,0]],masterActuator=None):
        ####################################################
        ### ACTUATOR

        # ANISOTROPY PARAMETER 
        # 0: YOUNG MODULUS TRANSVERSE
        # 1: YOUNG MODULUS LONGITUDINAL
        # 2: YAW
        # 3: ROLL

        ratio = [0.01,0.01,0.005,0.005] #[0.005,0.005,0.005,0.005] #[0.0001,0.0001,0.008,0.008] #[0.0001,0.0001,0.01,0.01] #[0.0001,0.0001,0.05,0.05]
        Ymax = 1000
        Ymin = 50
        # Ymin = 5
        angleMin = radians(-1000)
        angleMax = radians(1000)
        angleMin = radians(0)
        angleMax = radians(360)

        boundaries = [[Ymin,Ymax],[Ymin,Ymax],[angleMin,angleMax],[angleMin,angleMax]]
        # boundaries = [[angleMin,angleMax],[angleMin,angleMax]]
        if (masterActuator):
            # YOUNG MODULUS TRANSVERSE
            self.model.addObject('AnisotropyControlPointActuator', template='Vec3', name=self.name.value+'_actuator1',anisotropyParameter=anisotropyParameter,#cpNumber=[0],
                    forceField=self.model.Elasticity.getLinkPath(), maxVariationRatio=ratio,boundaries=boundaries,printLog=True,
                    masterActuator=masterActuator+"/master_actuator1")#,QPobjective=rootNode.QPInverseProblemSolver.objective) # 0.02
        else:
            # YOUNG MODULUS TRANSVERSE
            self.model.addObject('AnisotropyControlPointActuator', template='Vec3', name='master_actuator1',anisotropyParameter=anisotropyParameter,#cpNumber=[0],
                    forceField=self.model.Elasticity.getLinkPath(), maxVariationRatio=ratio, boundaries=boundaries,
                    targetParamValue=targetParamValue,#QPobjective=rootNode.QPInverseProblemSolver.getLinkPath()+".objective",
                    printLog=True) # 0.02
        return self

    def addController():
        self.model.addObject(Controller(name="c1_controller",goal=goal.goalMO))
        return self

    self.addModel = addModel
    self.attach = attach
    self.addEffector = addEffector
    self.addConstraint = addConstraint
    self.addGoal = addGoal
    self.addController = addController 

    return self