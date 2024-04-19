import Sofa
from splib.numerics import RigidDof, Quat, to_radians
from splib.animation import animate
from splib.constants import Key
from splib.numerics import Quat
import math
from math import *
import numpy as np

def quaternion_multiply(quaternion1, quaternion0):
    w0, x0, y0, z0 = quaternion0
    w1, x1, y1, z1 = quaternion1
    return np.array([-x1 * x0 - y1 * y0 - z1 * z0 + w1 * w0,
                     x1 * w0 + y1 * z0 - z1 * y0 + w1 * x0,
                     -x1 * z0 + y1 * w0 + z1 * x0 + w1 * y0,
                     x1 * y0 - y1 * x0 + z1 * w0 + w1 * z0], dtype=np.float64)

def setupanimation(actuators, step, angularstep, factor):
    """This functions is called repeatidely in an animation.
       It moves the actuators by translating & rotating them according to the factor
       value.
    """
    for actuator in actuators:
        rigid = RigidDof(actuator.dofs)
        rigid.translate(rigid.forward * step * factor)
        actuator.ServoMotor.angle += angularstep * factor

def rotateGoal(rigid, step, factor):
    if factor < 0.5:
        rigid.rotateAround([0,1,0],step)
    else:
        rigid.rotateAround([0,1,0],-step)

def euler_to_quaternion(roll, pitch, yaw):

    qx = np.sin(roll/2) * np.cos(pitch/2) * np.cos(yaw/2) - np.cos(roll/2) * np.sin(pitch/2) * np.sin(yaw/2)
    qy = np.cos(roll/2) * np.sin(pitch/2) * np.cos(yaw/2) + np.sin(roll/2) * np.cos(pitch/2) * np.sin(yaw/2)
    qz = np.cos(roll/2) * np.cos(pitch/2) * np.sin(yaw/2) - np.sin(roll/2) * np.sin(pitch/2) * np.cos(yaw/2)
    qw = np.cos(roll/2) * np.cos(pitch/2) * np.cos(yaw/2) + np.sin(roll/2) * np.sin(pitch/2) * np.sin(yaw/2)

    return [qx, qy, qz, qw]

class GoalController(Sofa.PythonScriptController):
    """This controller moves the goal position when the inverse control is activated
    """

    def __init__(self, goalNode,ellipsoid_pos,rigid_center_pos,workspaceVisu):
        self.name = "GoalController"
        self.activated = False
        self.time = 0
        self.dy = 0.1
        self.mo = goalNode.dofs
        self.rigidDOF = RigidDof(goalNode.dofs)
        self.stepsize = 0.5
        self.stepAngle = 5
        self.rotationAngle = to_radians(0.3) #0.5
        self.ellipsoid_trajectory = False

        self.ellipsoid_pos = ellipsoid_pos.dofs
        self.phantom = ellipsoid_pos

        self.rigid_center_pos = rigid_center_pos.dofs
        self.previousPos = [0,0,0,0,0,0,0]
        self.step = 6
        self.current_step = 1
        self.current_point_state = [True,True,True,True]
        self.workspaceVisu = workspaceVisu
        self.pos_count = 0
        self.rotate = False
        self.rotate_count = 0

    # def rotateOnItself(self,p,angle=[0,0,0]):
    #     q = Quat(p.value[0][3:])
    #     q_new = Quat.product(q.getInverse(), q)
    #     print([ i*180/math.pi for i in q.getInverse().getEulerAngles()])
    #     print("----------------")
    #     print(angle)
    #     # q_new.rotateFromEuler(angle)

    #     # Create a rotation object from Euler angles specifying axes of rotation
    #     rot_quat = Quat.createFromEuler(angle)
    #     q_new = Quat.product(q_new, rot_quat)
    #     print("rot_quat  ",rot_quat)

    #     print([ i*180/math.pi for i in q.getInverse().getEulerAngles()])
    #     # q = Quat.product(q,q_new)
    #     q = Quat([ i for i in quaternion_multiply(q,q_new)])
    #     print(quaternion_multiply(q,q_new))
    #     print([ i*180/math.pi for i in q.getInverse().getEulerAngles()])

    #     p.value = [p.value[0][0],p.value[0][1],p.value[0][2],q[0],q[1],q[2],q[3]]
    #     print("----------------")

    def rotateOnItself(self,p,angle=[0,0,0]):
        q = Quat(p.value[0][3:])
        q_new = Quat.product(q.getInverse(), q)
        q_new.rotateFromEuler(angle)
        q = Quat.product(q,q_new)
        p.value = [p.value[0][0],p.value[0][1],p.value[0][2],q[0],q[1],q[2],q[3]]

    def onKeyPressed(self, key):
        if key == Key.I:
            self.activated = True

        if key == Key.Z or key == Key.E:   
            if key == Key.Z:
                animate(rotateGoal, {"rigid": self.rigidDOF, "step": -self.rotationAngle}, duration=2) #,mode="loop")
            if key == Key.E:
                animate(rotateGoal, {"rigid": self.rigidDOF, "step": self.rotationAngle}, duration=2)#,mode="loop")
                # self.rigidDOF.rotateAround([1,0,0],to_radians(self.stepAngle))
                # self.rigidDOF.rotateAround([0,1,0],-self.rotationAngle)

        if key == Key.D:
            rigid_center = self.rigid_center_pos.position[0]
            angle = ToEulerAngles(rigid_center[3],rigid_center[4],rigid_center[5],rigid_center[6])
            print("-----------> ANGLE CENTER SIMULATED  : "+str(to_euler(angle[0]))+" | "+str(to_euler(angle[1]))+" | "+str(to_euler(angle[2])))



        p = self.rigidDOF.rigidobject.getData("position")
        if key == Key.T:
            self.rotateOnItself(p,[to_radians(self.stepAngle),0.,0.])
        if key == Key.G:
            self.rotateOnItself(p,[to_radians(-self.stepAngle),0.,0.])
        if key == Key.H:
            self.rotateOnItself(p,[0,to_radians(self.stepAngle),0.])
        if key == Key.Y:
            self.rotateOnItself(p,[0,to_radians(-self.stepAngle),0.])
        if key == Key.U:
            self.rotateOnItself(p,[0,0,to_radians(self.stepAngle)])
        if key == Key.J:
            self.rotateOnItself(p,[0.,0.,to_radians(-self.stepAngle)])
        if key == Key.K:
            self.ellipsoid_trajectory = not self.ellipsoid_trajectory
        if key == Key.L:
            # print('------------------------------------------------')
            # print(self.phantom.MeshTopology.position)
            for visu in self.workspaceVisu:
                print('------------------------------------------------')
                print(visu.position)
        else:
            pos = self.mo.position[0]

            if key == Key.uparrow:
                pos[1] += self.stepsize
            elif key == Key.downarrow:
                pos[1] -= self.stepsize

            if key == Key.leftarrow:
                pos[0] += self.stepsize
            elif key == Key.rightarrow:
                pos[0] -= self.stepsize

            if key == Key.plus:
                pos[2] += self.stepsize
            elif key == Key.minus:
                pos[2] -= self.stepsize

            self.mo.position = pos

    def rotateGoal(self,rotate_count):
        p = self.rigidDOF.rigidobject.getData("position")

        goalpos = [ i*180/math.pi for i in Quat(self.mo.position[0][3:]).getInverse().getEulerAngles()] 
        centerpos = [ i*180/math.pi for i in Quat(self.rigid_center_pos.position[0][3:]).getInverse().getEulerAngles()]           
        angle_offset = 3
        pos_offset = 1

        print(Quat(self.mo.position[0][3:]).getInverse(),Quat(self.rigid_center_pos.position[0][3:]).getInverse())
        print(goalpos,centerpos)

        print("rotation "+str(self.rotate_count))
        if abs(self.rigid_center_pos.position[0][0] - self.mo.position[0][0]) > pos_offset or abs(self.previousPos[0] - self.rigid_center_pos.position[0][0]) > pos_offset:
            if self.rotate_count <=2: self.current_point_state[1] = False
            # elif self.rotate_count <=4: self.current_point_state[2] = False
            # elif self.rotate_count <=6: self.current_point_state[3] = False
            print("Fail rotation "+str(self.rotate_count))

        if abs(self.rigid_center_pos.position[0][1] - self.mo.position[0][1]) > pos_offset or abs(self.previousPos[1] - self.rigid_center_pos.position[0][1]) > pos_offset:   
            if self.rotate_count <=2: self.current_point_state[1] = False
            # elif self.rotate_count <=4: self.current_point_state[2] = False
            # elif self.rotate_count <=6: self.current_point_state[3] = False
            print("Fail rotation "+str(self.rotate_count))

        if abs(self.rigid_center_pos.position[0][2] - self.mo.position[0][2]) > pos_offset or abs(self.previousPos[2] - self.rigid_center_pos.position[0][2]) > pos_offset:
            if self.rotate_count <=2: self.current_point_state[1] = False
            # elif self.rotate_count <=4: self.current_point_state[2] = False
            # elif self.rotate_count <=6: self.current_point_state[3] = False
            print("Fail rotation "+str(self.rotate_count))

        if abs(centerpos[0]-goalpos[0]) > angle_offset or abs(centerpos[1]-goalpos[1]) > angle_offset or abs(centerpos[2]-goalpos[2]) > angle_offset:
            print("Fail rotation due to Translation "+str(self.rotate_count))

            if self.rotate_count <=2: self.current_point_state[1] = False
            # if self.rotate_count <=4: self.current_point_state[2] = False
            # if self.rotate_count <=6: self.current_point_state[3] = False

        if rotate_count == 0:
            # print("rot0")
            # print(goalpos,centerpos)
            # self.rotateOnItself(p,[to_radians(self.stepAngle),0.,0.])
            self.rotateOnItself(p,[0.,to_radians(self.stepAngle),0.])

            # goalpos = [ i*180/math.pi for i in Quat(self.mo.position[0][3:]).getInverse().getEulerAngles()] 
            # centerpos = [ i*180/math.pi for i in Quat(self.rigid_center_pos.position[0][3:]).getInverse().getEulerAngles()]

            # print(Quat(self.mo.position[0][3:]).getInverse(),Quat(self.rigid_center_pos.position[0][3:]).getInverse())
            # print(goalpos,centerpos)
            return

            # print(self.stepAngle)
        if rotate_count == 1:
            # self.rotateOnItself(p,[to_radians(-2*self.stepAngle),0.,0.])
            self.rotateOnItself(p,[0.,to_radians(-2*self.stepAngle),0.])

            return
            # print(-2*self.stepAngle)
        if rotate_count == 2:
            # self.rotateOnItself(p,[to_radians(self.stepAngle),0.,0.])
            self.rotateOnItself(p,[0.,to_radians(self.stepAngle),0.])

            self.rotate_count = 0
            self.rotate = False
            return
        #     self.rotateOnItself(p,[to_radians(self.stepAngle),0.,0.])
        #     self.rotateOnItself(p,[0,to_radians(self.stepAngle),0.])
        # if rotate_count == 3:
        #     self.rotateOnItself(p,[0,to_radians(-2*self.stepAngle),0.])
        # if rotate_count == 4:
        #     self.rotateOnItself(p,[0,to_radians(self.stepAngle),0.])
        #     self.rotateOnItself(p,[0,0,to_radians(self.stepAngle)])
        # if rotate_count == 5:
        #     self.rotateOnItself(p,[0.,0.,to_radians(-2*self.stepAngle)])
        # if rotate_count == 6:
        #     self.rotateOnItself(p,[0,0,to_radians(self.stepAngle)])
        #     self.rotate_count = 0
        #     self.rotate = False

    def onBeginAnimationStep(self, dt):
        if self.ellipsoid_trajectory:
            # print(self.current_step,self.step)
            self.current_step += 1
            if self.rotate:
                if self.current_step % self.step == self.step-2 and self.pos_count >= 1:
                    self.previousPos = self.rigid_center_pos.position[0]
                if self.rotate_count == 0 and self.current_step % self.step == self.step-1:
                    offset = 1
                    self.current_point_state[0] = False
                    if abs(self.rigid_center_pos.position[0][0] - self.mo.position[0][0]) < offset and abs(self.previousPos[0] - self.rigid_center_pos.position[0][0]) < offset:
                        if abs(self.rigid_center_pos.position[0][1] - self.mo.position[0][1]) < offset and abs(self.previousPos[1] - self.rigid_center_pos.position[0][1]) < offset:   
                            if abs(self.rigid_center_pos.position[0][2] - self.mo.position[0][2]) < offset and abs(self.previousPos[2] - self.rigid_center_pos.position[0][2]) < offset:
                                self.current_point_state[0] = True
                                # print("TRUE  : "+str(tmp))
                            else:
                                print("Translation ERROR")
                        else:
                            print("Translation ERROR")
                    else:
                        print("Translation ERROR")

                # goalpos = [ i*180/math.pi for i in Quat(self.mo.position[0][3:]).getInverse().getEulerAngles()] 
                # centerpos = [ i*180/math.pi for i in Quat(self.rigid_center_pos.position[0][3:]).getInverse().getEulerAngles()]           
                # print("befrore   ",goalpos,centerpos)
                if self.current_step % self.step == 0:
                    self.rotateGoal(self.rotate_count)
                    self.rotate_count += 1

            else:
                if self.current_step % self.step == 0:
                    if self.pos_count >= 1:
                        tmp_visu = -1
                        if all( i is True for i in self.current_point_state ):                                        tmp_visu = 0 # DOF == 3T + 3R
                        elif self.current_point_state[0] and self.current_point_state[1:].count(False) == 1  :        tmp_visu = 1 # DOF == 3T + 2R
                        elif self.current_point_state[0] and self.current_point_state[1:].count(False) == 2 :         tmp_visu = 2 # DOF == 3T + 1R 
                        elif self.current_point_state[0] and all( i is False for i in self.current_point_state[1:] ): tmp_visu = 3 # DOF == 3T 
                        else: tmp_visu = 4
                        print("-----------------------------")
                        print(self.current_point_state,tmp_visu)
                        if tmp_visu != -1:
                            print('VISU')
                            tmp = self.workspaceVisu[tmp_visu].position
                            tmp.append(self.ellipsoid_pos.position[self.pos_count-1])
                            self.workspaceVisu[tmp_visu].position = tmp
                            self.workspaceVisu[tmp_visu].init()

                        self.current_point_state = [True,True,True,True]

                    # print(self.pos_count)
                    # print("-----------------------------")
                    # print(self.mo.position)
                    # print(self.ellipsoid_pos.position[self.pos_count])
                    # tmp = self.mo.position[0]
                    # tmp[:3] = self.ellipsoid_pos.position[self.pos_count]
                    # print(tmp)
                    tmp = self.mo.position[0]
                    tmp[:3] = self.ellipsoid_pos.position[self.pos_count]  
                    self.mo.position = tmp
                    self.pos_count += 1
                    self.rotate_count = 0 
                    self.rotate = True
                # print(self.current_step % self.step,self.pos_count,self.step)
                # if self.current_step % self.step == self.step-2 and self.pos_count >= 1:
                #     self.previousPos = self.rigid_center_pos.position[0]
                # if self.current_step % self.step == self.step-1 and self.pos_count >= 1:
                #     offset = 1
                #     if not self.rotate: 
                #         self.current_point_state[0] = False
                #         if abs(self.rigid_center_pos.position[0][0] - self.mo.position[0][0]) < offset and abs(self.previousPos[0] - self.rigid_center_pos.position[0][0]) < offset:
                #             if abs(self.rigid_center_pos.position[0][1] - self.mo.position[0][1]) < offset and abs(self.previousPos[1] - self.rigid_center_pos.position[0][1]) < offset:   
                #                 if abs(self.rigid_center_pos.position[0][2] - self.mo.position[0][2]) < offset and abs(self.previousPos[2] - self.rigid_center_pos.position[0][2]) < offset:
                #                     self.current_point_state[0] = True
                #                     # print("TRUE  : "+str(tmp))
                #                 else:
                #                     print("Translation ERROR")
                #             else:
                #                 print("Translation ERROR")
                #         else:
                #             print("Translation ERROR")
    #     if self.activated:
    #         self.time = self.time+dt

    #     if self.time >= 1:
    #         self.time = 0;
    #         self.dy = -self.dy

    #     pos = [self.mo.position[0][0], self.mo.position[0][1], self.mo.position[0][2]]
    #     pos[1] += self.dy
    #     self.mo.position = [[pos[0], pos[1], pos[2], 0, 0, 0, 1]]


class DirectController(Sofa.PythonScriptController):
    """This controller has two role:
       - if user press up/left/right/down/plus/minus the servo motor angle
         is changed.
       - if user press A an animation is started to move the motor to the physical position
         they are occupying in the real robot.
    """

    def __init__(self, node, actuators, serialportctrl):
        self.name = "TripodController"
        self.stepsize = 0.1
        self.actuators = actuators
        self.serialportctrl = serialportctrl

    def onKeyPressed(self, key):
        if key == Key.A and self.serialportctrl.state == "init":
            self.serialportctrl.state = "no-comm"
            animate(setupanimation, {"actuators": self.actuators, "step": 3.0, "angularstep": -0.14}, duration=0.2) # -0.14

        # Inclusion of the keystroke to start data sending = establishing communication ('comm')
        if key == Key.B and self.serialportctrl.state == "no-comm":
            self.serialportctrl.state = "comm"


class InverseController(Sofa.PythonScriptController):
    """This controller has two role:
       - if user press up/left/right/down/plus/minus the servo motor angle
         is changed.
       - if user press A an animation is started to move the motor to the physical position
         they are occupying in the real robot.
    """

    def __init__(self, node, nodeGoal, nodeEffector, nodeActuators, nodeDofRigid, nodeTripod, serialport=None, servomotors=None,limitLow=None,limitHigh=None):
        self.name = "InverseController"
        self.node = node
        self.nodeGoal = nodeGoal
        self.nodeEffector = nodeEffector
        self.nodeActuators = nodeActuators
        self.nodeDofRigid = nodeDofRigid
        self.nodeTripod = nodeTripod
        self.serialport = serialport
        # self.serialport.packetOut = [150, 150, 150]
        self.state = "init"
        self.actuators = servomotors
        self.activate = False
        self.limitHigh = limitHigh
        self.limitLow = limitLow 
    def onKeyPressed(self, key):
        if key == Key.I:
            self.activate = 1

    def onBeginAnimationStep(self,dt):
        self.nodeActuators.activated = bool(self.activate)
        self.node.Simulation.rigid_center.posEffector.activated = bool(self.activate)
        self.nodeActuators.init()
        self.node.Simulation.rigid_center.posEffector.init()

    def onEndAnimationStep(self,dt):

        # if self.state == "init":
        #     return

        # if self.state == "no-comm":
        #     return

        if(self.nodeActuators.activated):
            W_R_Dof = [[0]*4,[0]*4,[0]*4];
            Angles = [0]*3;

            for i in range(0,4):
                W_R_Dof[0][i] = self.nodeDofRigid.dofs.position[0][i+3]
                W_R_Dof[1][i] = self.nodeDofRigid.dofs.position[1][i+3]
                W_R_Dof[2][i] = self.nodeDofRigid.dofs.position[2][i+3]

            # ActuatedArm0
            Angles[0] = Quat(Quat(W_R_Dof[0]).getInverse()).getEulerAngles()[0]*180/math.pi
            # if Angles[0] < 0:
            armRotation = 2.0943951
            # else:
            #     armRotation = 2.0943951
            # ActuatedArm1
            q = Quat(Quat(W_R_Dof[1]).getInverse())
            q.rotateFromEuler([0.,armRotation,0.]) # 120 degree
            Angles[1] = q.getEulerAngles()[0]*180/math.pi
            # ActuatedArm2
            q = Quat(Quat(W_R_Dof[2]).getInverse())
            q.rotateFromEuler([0.,armRotation*2,0.]) # # 240 degree
            Angles[2] = q.getEulerAngles()[0]*180/math.pi
            # print(Angles) 
            Angles = [abs(int(angle))+self.limitLow[i] for i, angle in enumerate(Angles)]

            for i in range(3):
                if Angles[i] < self.limitLow[i]:
                    Angles[i] = self.limitLow[i]
                    # print("----------------------->     "+str(self.limitLow[i]))
                if Angles[i] > self.limitHigh[i]: 
                    Angles[i] = self.limitHigh[i]
                    # print("----------------------->     "+str(self.limitHigh[i]))

                Angles[i] = self.limitHigh[i] + self.limitLow[i] - Angles[i]

            # print(Angles)
            # print("===================")
            # The controller board of the real robot receives `AnglesDeg` values
            if(self.serialport):
                self.serialport.packetOut = Angles

def ToEulerAngles(x,y,z,w):

    angles = []

    # roll (x-axis rotation)
    sinr_cosp = +2.0 * (w * x + y * z)
    cosr_cosp = +1.0 - 2.0 * (x * x + y * y)
    angles.append(atan2(sinr_cosp, cosr_cosp))

    # pitch (y-axis rotation)
    sinp = +2.0 * (w * y - z * x)
    if (fabs(sinp) >= 1):
        angles.append(copysign(M_PI / 2, sinp)) # use 90 degrees if out of range
    else:
        angles.append(asin(sinp))

    # yaw (z-axis rotation)
    siny_cosp = +2.0 * (w * z + x * y)
    cosy_cosp = +1.0 - 2.0 * (y * y + z * z); 
    angles.append(atan2(siny_cosp, cosy_cosp))

    return angles

def to_euler(angle):

    return angle * 180/math.pi

class CenterController(Sofa.PythonScriptController): #Sofa.Core.Controller):

    def __init__(self, node):
        self.node = node

    def onEndAnimationStep(self, dt):

        rigid_center = self.node.dofs.position[0]
        angle = ToEulerAngles(rigid_center[3],rigid_center[4],rigid_center[5],rigid_center[6])
        print("----------->  "+str(to_euler(angle[0])-60)+" | "+str(to_euler(angle[1]))+" | "+str(to_euler(angle[2])+90))
