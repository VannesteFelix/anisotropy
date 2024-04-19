import os, sys

from splib.numerics import sin, cos, to_radians
from stlib.physics.deformable import ElasticMaterialObject
from stlib.physics.collision import CollisionMesh
from splib.objectmodel import SofaPrefab, SofaObject
from stlib.components import addOrientedBoxRoi
from splib.numerics import vec3
from splib.numerics.quat import Quat
from splib.numerics import *
from stlib.solver import DefaultSolver
from stlib.scene import Node

path = os.path.dirname(os.path.abspath(__file__))
path_utility = path+"/../../utility/"
sys.path.insert(0,path_utility)
from rigidification import Rigidify
from highOrderShape import HighOrderShape

from actuatedarm import ActuatedArm # from . actuatedarm import ActuatedArm
from tutorial import * # from . tutorial import *

def to_euler(angle):

    return angle * 180/math.pi

def ElasticBody(parent,name="ElasticMaterialObject",ho=False,**kwargs):

    body = parent.createChild("ElasticBody")

    if ho:
        e = HighOrderShape(body,name=name,**kwargs)
        e = e.highOrder_node
    else:
        e = ElasticMaterialObject(body,name=name,**kwargs)
        # e.forcefield.drawHeterogeneousTetra = True
        # e.forcefield.drawAsEdges = True

    # visual = body.createChild("Visual")
    # visual.createObject("MeshSTLLoader", name="loader", filename="mesh/tripod_mid.stl")
    # visual.createObject("OglModel", name="renderer", src="@loader", color=[1.0, 1.0, 1.0, 0.5],
    #                     rotation=[90, 0, 0], translation=[0, 30, 0])

    # visual.createObject("BarycentricMapping",
    #                     input=e.dofs.getLinkPath(),
    #                     output=visual.renderer.getLinkPath())
    
    return body


@SofaPrefab
class Tripod(SofaObject):

    def __init__(self, parent, name="Tripod", radius=66, numMotors=3, angleShift=180.0,ho=False,direct=False,inverse=False,**kwargs):
        self.direct = direct
        self.inverse = inverse
        self.node = parent.createChild(name)
        ElasticBody(self.node,ho=ho,**kwargs)
        
        self.translation = kwargs['translation']
        self.rotation = kwargs['rotation']

        self.ho = ho
        self.body = self.node.ElasticBody.ElasticMaterialObject

        dist = radius
        numstep = numMotors
        self.actuatedarms = []
        self.trsform = []
        for i in range(0, numstep):
            name = "ActuatedArm"+str(i)
            tr, eulerRotation = self.__getTransform(i, numstep, angleShift, radius, dist-0.8)
            #print(tr,eulerRotation)
            arm = ActuatedArm(self.node, name=name,
                              translation=tr, eulerRotation=eulerRotation)
            self.actuatedarms.append(arm)
            # Add limits to angle that correspond to limits on real robot
            arm.ServoMotor.minAngle = -2.0225
            arm.ServoMotor.maxAngle = -0.0255
            if i == 0:
                self.trsform.append(self.node.ActuatedArm0.ServoMotor.BaseFrame.dofs.position)
            else:
        
                angle = i*360/3
                q = Quat(self.node.ActuatedArm0.ServoMotor.BaseFrame.dofs.position[0][3:])
                q.rotateFromEuler([0,to_radians(angle),0])
                pos = self.node.ActuatedArm0.ServoMotor.BaseFrame.dofs.position
        
                pos[0][:3] = arm.ServoMotor.BaseFrame.dofs.position[0][:3]
                pos[1][:3] = arm.ServoMotor.BaseFrame.dofs.position[0][:3]
                pos[0][3:] = q
                pos[1][3:] = q
                self.trsform.append(pos)    
        print(self.trsform)
        self.node.ActuatedArm1.ServoMotor.BaseFrame.dofs.position = self.trsform[1]
        self.node.ActuatedArm2.ServoMotor.BaseFrame.dofs.position = self.trsform[2]

        self.__attachToActuatedArms(radius, numMotors, angleShift)
        if inverse:
            self.__addActuation(numMotors,*inverse)
        
        
        #for arm in self.actuatedarms:
        #    print(arm.node)
        #    arm.node.ServoArm.addChild(self.node.RigidifiedStructure.RigidParts)
        #self.node.RigidifiedStructure.RigidParts.addChild(arm.node)

    def __getTransform(self, index, numstep, angleShift, radius, dist,tr=0):
        fi = float(index)
        fnumstep = float(numstep)
        angle = fi*360/fnumstep
        angle2 = fi*360/fnumstep+angleShift
        eulerRotation = [0, angle, 0]
        translation = [dist*sin(to_radians(angle2)), -1.35+tr, dist*cos(to_radians(angle2))]

        eulerRotation = [i+j for i,j in zip(eulerRotation,self.rotation)]
        translation = transformPositions([translation], translation=self.translation, eulerRotation=self.rotation, scale=[1.0,1.0,1.0])[0] #self.rotation

        return translation, self.rotation

    def addCollision(self):
        # CollisionMesh(self.body, surfaceMeshFileName="mesh/collision_arm.stl", name="CollisionModel", translation=[0.0, 28, 0.0], rotation=[90, 180, 0], collisionGroup=1)
        # CollisionMesh(self.body, surfaceMeshFileName="mesh/collision_arm.stl", name="CollisionModel", translation=[0.0, 28, 0.0], rotation=[90, 300, 0], collisionGroup=1)
        CollisionMesh(self.body, surfaceMeshFileName="mesh/collision_arm_little.stl", name="CollisionModel", translation=[0.0, 28, 0.0], rotation=[90, 240, 0], collisionGroup=1)
 
        # CollisionMesh(self.body, surfaceMeshFileName="mesh/tripod_low.stl", name="CollisionModel", translation=[0.0, 28, 0.0], rotation=[90, 0, 0], collisionGroup=1)
        # CollisionMesh(self.body, surfaceMeshFileName="mesh/tripod_low.stl", name="CollisionModel", translation=[0.0, 28, 0.0], rotation=[90, 0, 0], collisionGroup=1)

        # CollisionMesh(self.body, surfaceMeshFileName="mesh/tripod_low.stl", name="CollisionModel", translation=[0.0, 28, 0.0], rotation=[90, 0, 0], collisionGroup=1)

        # CollisionMesh(self.body, surfaceMeshFileName="mesh/tripod_whole_collision.stl", name="CollisionModel", translation=[0.0, 28, 0.0], rotation=[90, 0, 0], collisionGroup=[1,2,3])
        # CollisionMesh(self.body, surfaceMeshFileName="mesh/tripod_arm_collision.stl", name="CollisionModel_0", translation=[0.0, 28, 0.0], rotation=[90, 0, 0], collisionGroup=1)
        # CollisionMesh(self.body, surfaceMeshFileName="mesh/tripod_arm_collision.stl", name="CollisionModel_1", translation=[0.0, 28, 0.0], rotation=[90, 120, 0], collisionGroup=2)
        # CollisionMesh(self.body, surfaceMeshFileName="mesh/tripod_arm_collision.stl", name="CollisionModel_2", translation=[0.0, 28, 0.0], rotation=[90, -120, 0], collisionGroup=3)

        for i,arm in enumerate(self.actuatedarms):
        #     CollisionMesh(arm.ServoMotor.ServoBody,
        #                   surfaceMeshFileName="mesh/servo_collision.stl",
        #                   name="TopServoCollision", mappingType='RigidMapping')
            if i == 2:
                CollisionMesh(arm.ServoArm,
                             surfaceMeshFileName="mesh/SG90_servoarm.stl",translation=[0.0, -25, 0.0],
                             name="CollisionMeshAuto2",mappingType='RigidMapping',collisionGroup=2)

    def __attachToActuatedArms(self, radius=66, numMotors=3, angleShift=180.0):
        deformableObject = self.body

        dist = radius
        numstep = numMotors
        groupIndices = []
        frames = []

        noInit=False        
        # Rigidify the deformable part at extremity to attach arms
        if self.ho:
            noInit=True

            deformableObject.loader.init()
            deformableObject.Container1.init()
            deformableObject.GeomAlgo.init()
            deformableObject.dofs.init()

        deformableObject = self.node.ElasticBody.ElasticMaterialObject

        for i in range(0, numstep):
            # translation, eulerRotation = self.__getTransform(i, numstep, angleShift, radius, dist-2)
            # box = addOrientedBoxRoi(self.node, position=deformableObject.dofs.getData("rest_position"), name="BoxROI"+str(i),
            #                         translation=vec3.vadd(translation, [0.0, 25.0, 0.0]),
            #                         eulerRotation=eulerRotation, scale=[45, 15, 30])
        
            # box.drawBoxes = True
            # box.init()
            # groupIndices.append([ind[0] for ind in box.indices])
        
        
            translation, eulerRotation = self.__getTransform(i, numstep, angleShift, radius, dist,25)
            frames.append(translation + self.trsform[i][0][3:])
        #print(frames)


       ## ADD CENTER
        # o = deformableObject.createObject("SphereROI", name="roi", template="Rigid3",
        #             position=deformableObject.dofs.getData("rest_position"),
        #             centers=[0.0, 28, 0.0], radii=[7.5], drawSphere=True)
        # o.init()
        # tmp = [i[0] for i in o.indices] 
        # print(tmp)
        translation = [0.0, 28, 0.0]
        translation = transformPositions([translation], translation=self.translation, eulerRotation=self.rotation, scale=[1.0,1.0,1.0])[0]
        frames.append( translation + self.trsform[0][0][3:])
        print(frames)
        ####################

        # print(groupIndices)
        groupIndices = [[38, 39, 42, 43, 205, 206, 207, 208, 209, 210, 211, 212, 221, 222, 237, 238, 239, 240, 241, 242, 243, 244, 253, 254, 258, 259, 712, 713, 719, 732, 733, 739, 740, 743, 748, 751, 768, 770, 777, 778, 779, 783, 784, 814, 815, 821, 834, 835, 841, 842, 845, 850, 853, 870, 872, 879, 880, 881, 885, 886, 932, 933, 940, 943, 944, 945, 946, 947, 948, 949, 950, 951, 952, 953, 961],
                        [23, 24, 32, 33, 103, 104, 105, 106, 107, 108, 109, 110, 119, 120, 165, 166, 167, 168, 169, 170, 171, 172, 181, 182, 255, 257, 386, 387, 393, 408, 409, 415, 416, 420, 426, 428, 445, 447, 454, 455, 456, 460, 461, 519, 520, 526, 541, 542, 548, 549, 553, 559, 561, 578, 580, 587, 588, 589, 593, 594, 897, 898, 899, 908, 909, 910, 911, 912, 913, 914, 915, 916, 917, 918, 929],
                        [5, 6, 26, 27, 59, 60, 61, 62, 63, 64, 65, 66, 75, 76, 131, 132, 133, 135, 136, 137, 138, 139, 140, 141, 150, 151, 274, 275, 281, 294, 295, 301, 302, 305, 310, 313, 330, 332, 339, 340, 341, 345, 346, 472, 473, 476, 483, 484, 485, 486, 487, 488, 489, 490, 491, 494, 495, 499, 617, 618, 624, 637, 638, 644, 645, 648, 653, 656, 673, 675, 682, 683, 684, 688, 689],
                        [1, 2, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 46, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 355, 356, 357, 358, 359, 360, 361, 362, 363, 364, 365, 366, 367, 368, 369, 370, 371, 372, 373, 374, 375, 376, 377, 378, 419, 424, 444, 537, 552, 557, 577, 966, 967]]

        
        if self.inverse:
            frames.append([0, 30, 0, 0, 0, 0])
            groupIndices.append([1, 2, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 46, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 355, 356, 357, 358, 359, 360, 361, 362, 363, 364, 365, 366, 367, 368, 369, 370, 371, 372, 373, 374, 375, 376, 377, 378])

            #center = self.node.createObject("SphereROI", name="roi", template="Rigid3",
            #                            position=deformableObject.dofs.getData("rest_position"),
            #                            centers=[0, 28, 0], radii=[6], drawSphere=True)
            #center.init()
            #groupIndices.append([ i[0] for i in center.indices])
            
        rigidifiedstruct = Rigidify(self.node, deformableObject, groupIndices=groupIndices, frames=frames, name="RigidifiedStructure",noInit=noInit,ho=self.ho)

        if not self.direct:    
            rigidifiedstruct.DeformableParts.createObject("UncoupledConstraintCorrection")
            rigidifiedstruct.RigidParts.createObject("UncoupledConstraintCorrection")
    
        # Use this to activate some rendering on the rigidified object ######################################
        setData(rigidifiedstruct.RigidParts.dofs, showObject=True, showObjectScale=10, drawMode=2)
        # setData(rigidifiedstruct.RigidParts.RigidifiedParticules.dofs, showObject=True, showObjectScale=1,drawMode=1, showColor=[1., 1., 0., 1.])
        # setData(rigidifiedstruct.DeformableParts.dofs, showObject=True, showObjectScale=1, drawEdges=1, showColor=[0., 0., 1., 1.])
        #####################################################################################################

        # Attach arms
        for i in range(0, numstep):
           rigidifiedstruct.RigidParts.createObject('RestShapeSpringsForceField', name="rssff"+str(i),
                                                    points=i,
                                                    external_rest_shape=self.actuatedarms[i].servoarm.dofs,
                                                    stiffness='1e16', angularStiffness='1e7')
        # print(self.actuatedarms[i].servoarm.dofs)
        # rigidifiedstruct.RigidParts.createObject('RestShapeSpringsForceField', name="rssff"+str(numstep),
        #                                  points=numstep,
        #                                  stiffness='1e12', angularStiffness='0')
        

        # NEW CONFIG ---------- NOT WORKING
        # for arm in self.actuatedarms:
        #     print(arm.node)
        #     arm.node.ServoArm.addChild(self.node.RigidifiedStructure.RigidParts)
        # for i in range(3):
        #     rigidifiedstruct.RigidParts.createObject('RigidRigidMapping',
        #                            name="mapping", input=self.actuatedarms[i].node.ServoArm.dofs, index=0)

    def __addActuation(self,numMotors=3,goalNode=None,use_orientation=True):
        
        actuators = self.RigidifiedStructure.RigidParts.createChild('actuators')
        print(self.actuatedarms[0].ServoMotor.Angle.getLinkPath())
        #actuators.createObject('SlidingActuator', name="SlidingActuator0", template='Rigid3d', direction='0 0 0 1 0 0' , indices=0, maxForce='150000', minForce='-150000') #,context=self.actuatedarms[0].ServoMotor.Angle.getLinkPath())
        #actuators.createObject('SlidingActuator', name="SlidingActuator1", template='Rigid3d', direction='0 0 0 '+str(cos(4*math.pi/3))+' 0 '+str(sin(4*math.pi/3)) , indices=1, showDirection='1', showVisuScale='100', maxForce='150000', minForce='-150000')
        #actuators.createObject('SlidingActuator', name="SlidingActuator2", template='Rigid3d', direction='0 0 0 '+str(cos(2*math.pi/3))+' 0 '+str(sin(2*math.pi/3)) , indices=2, showDirection='1',  showVisuScale='100', maxForce='150000', minForce='-150000')
    
        if goalNode==None:
            actuators.createObject('PositionEffector', name='effector', template='Rigid3d', useDirections='1 1 1 0 0 0', indices='3', effectorGoal="10 40 0", limitShiftToTarget=True, maxShiftToTarget=5 )
        elif use_orientation:
            # Effector = actuators.createObject('PositionEffector', name='effector', template='Rigid3d', useDirections='0 1 0 1 0 1', indices='3', effectorGoal=goalNode.goalMO.getLinkPath()+".position", limitShiftToTarget=True, maxShiftToTarget=5)
            actuators.createObject('PositionEffector', name='effector', template='Rigid3d', useDirections='0 0 0 1 1 1', indices='3', effectorGoal=goalNode.goalMO.getLinkPath()+".position", limitShiftToTarget=True, maxShiftToTarget=5)
        else:
            actuators.createObject('PositionEffector', name='effector', template='Rigid3d', useDirections='1 1 1 0 0 0', indices='3', effectorGoal=goalNode.goalMO.getLinkPath()+".position", limitShiftToTarget=True, maxShiftToTarget=5)
    
        actuators.activated = 0
