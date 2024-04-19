#Sofa import
import Sofa
#Splib import
from splib3.numerics import RigidDof, to_radians
from splib3.animation import animate
from splib3.constants import Key
#Local import
from tutorial import *
#Other import
import math



def setupanimation(actuators, step, angularstep, factor):
    """This function is called repeatidely in an animation.
       It moves the actuators by translating & rotating them according to the factor
       value.
    """
    for actuator in actuators:
        rigid = RigidDof( actuator.ServoMotor.BaseFrame.dofs )
        rigid.setPosition( rigid.rest_position + rigid.forward * (step+5) * factor ) # (step+5)
        actuator.angleIn = angularstep * factor
        

def moveMotorSynchronously(actuators, angularstart,angularend,offset,down, factor):
    """This function is called repeatidely in an animation.
       It moves the actuators by translating & rotating them according to the factor
       value.
    """
    angularstep = -1 * (angularstart - angularend)
    for actuator in actuators:
        actuator.angleIn = angularstart + angularstep * factor
    if offset != None:
        offset.calibrate = False #True
        if down:
            test = [0.0, 0.0, 0.0] #[-0.3, -0.2, -0.2]
        else:
            test = [0.0, 0.0, 0.0]
        for i,actuator in enumerate(actuators):
            angularstep = -1 * (angularstart - (angularend+test[i]))
            offset.manual[i] = angularstart + angularstep * factor
            # print("calibrating "+str(offset.manual[i]))
        # if factor == 1:
        #     print("succes!!!")
        #     offset.calibrate = False
    else:
        test = [0.0, 0.0, 0.0]
        for i,actuator in enumerate(actuators):
            angularstep = -1 * (angularstart - (angularend+test[i]))
            # offset.manual[i] = angularstart + angularstep * factor

def saveTripodPosition(actuators, step, angularstep, factor):
        t = []
        for actuator in actuators:
                t.append(actuator.ServoMotor.BaseFrame.dofs.getData("position"))
                t.append(actuator.ServoMotor.getData("angleIn"))

        dumpPosition(t, "tripodRestPosition.json")

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


from splib.numerics.vec3 import Vec3
import math
def spiral(container,mecha):

    directions = []
    
    for tetrahedra in container.tetrahedra:
        barycenter = Vec3()
        # print("TETRAHEDRA : "+str(tetrahedra))
        for pointIndice in tetrahedra:
            # print("TETRAHEDRA INDICES : "+str(pointIndice))
            # print("MECHA POSITION : "+str(mecha.position[pointIndice]))

            barycenter += mecha.position[pointIndice]

        barycenter/=4

        direction = (math.atan2(barycenter[0],barycenter[2])*180/math.pi) #-90 #+ 50 #(150.0 +)/360
        if direction < 0:
            direction += 360 
        if direction > 360:
            direction -= 360
        directions.append(direction)
    
    return directions

def rotateGoal(rigid, step, factor):
    if factor < 0.5:
        rigid.rotateAround([0,0,1],step)
    else:
        rigid.rotateAround([0,0,1],-step)

def RepresentsInt(s):
    try: 
        int(s)
        return True
    except ValueError:
        return False

class GoalController(Sofa.PythonScriptController):
    """This controller moves the goal position when the inverse control is activated
    """

    def __init__(self, goalNode,nodeActuators):
        self.name = "GoalController"
        self.activated = False
        self.time = 0
        self.dy = 0.1
        self.mo = goalNode.goalMO
        self.rigidDOF = RigidDof(goalNode.goalMO)
        self.stepsize = 0.5
        self.rotationAngle = to_radians(0.5)
        self.nodeActuators = nodeActuators


    def onKeyPressed(self, key):
        if key == Key.I:
            self.activated = True
            for actuator in self.nodeActuators:
                actuator.activated = bool(self.activated)
                actuator.init()
            #self.nodeActuators.activated = bool(self.activated)
            #self.nodeActuators.init()
            print("##########################################")
            print("-----------> inverse control activated")
            print("##########################################")

        if self.activated :
            if key == Key.J or key == Key.K:   
                if key == Key.J:
                    print("here")
                    animate(rotateGoal, {"rigid": self.rigidDOF, "step": self.rotationAngle}, duration=2,mode="loop")
                if key == Key.K:
                    self.rigidDOF.rotateAround([0,1,0],to_radians(-14))
                    # self.rigidDOF.rotateAround([0,1,0],-self.rotationAngle)
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

class DirectController(Sofa.PythonScriptController):
    """This controller has two role:
       - if user press up/left/right/down/plus/minus the servo motor angle
         is changed.
       - if user press A an animation is started to move the motor to the physical position
         they are occupying in the real robot.
    """

    def __init__(self, node, actuators):
        self.name = "TripodController"
        self.stepsize = 0.01
        self.actuators = actuators
        self.angleInit = -1.2  #-1.2 #-1.4965
        self.angleMiddle = -0.87  
        self.angleTop =  -0.35 #-0.54 #-0.6 #-0.35 #-0.6 #-0.4 # -0.7 #-0.85 # -0.6


    def onKeyPressed(self, key):
        self.initTripod(key)

    def initTripod(self, key):
        if key == Key.A:
            animate(setupanimation, 
                    {"actuators": self.actuators, "step": 35.0, "angularstep": self.angleInit},
                    duration=2,
                    onDone=saveTripodPosition)

class CenterController(Sofa.PythonScriptController): #Sofa.Core.Controller):

    def __init__(self, node):
        self.node = node

    def onEndAnimationStep(self, dt):

        rigid_center = self.node.dofs.position[0]
        angle = ToEulerAngles(rigid_center[3],rigid_center[4],rigid_center[5],rigid_center[6])
        print("----------->  "+str(to_euler(angle[0]))+" | "+str(to_euler(angle[1]))+" | "+str(to_euler(angle[2])))

class MultipleTripodController(Sofa.PythonScriptController): 

    def __init__(self, node, tripods):
        self.node = node
        self.tripods = tripods

        self.currentControlledTripod = 0
        self.currentActuator = [tr.node for tr in self.tripods[0].actuatedarms]


        self.angleInit = -1.2  #-1.2 #-1.4965
        self.angleMiddle = -0.87  
        self.angleTop =  -0.35 #-0.54 #-0.6 #-0.35 #-0.6 #-0.4 # -0.7 #-0.85 # -0.6
        self.stepsize = 0.01


    def onKeyPressed(self, key):
        if key == Key.A:
            self.initTripod(key)
        self.animateTripod(key)

    def initTripod(self, key):
        for tripod in self.tripods:
            animate(setupanimation, 
                    {"actuators": [tr.node for tr in tripod.actuatedarms], "step": 35.0, "angularstep": self.angleInit},
                    duration=1,
                    onDone=saveTripodPosition)

    def moveMotors(self,keyToActivate,angularstart,angularend,key,actuators,offset=None,down=False):
        if key == keyToActivate:
            animate(moveMotorSynchronously, 
                    {"actuators": actuators,"angularstart": angularstart, "angularend": angularend,'offset':offset, 'down':down}, 
                    duration=1)

    def animateTripod(self, key):

        # self.moveMotors(Key.E,self.angleTop,self.angleInit,key)
        if RepresentsInt(key):
            if int(key) in range(10):
                print("New current controlled Tripod --------->  "+key)
                self.currentControlledTripod = int(key)
                self.currentActuator = [tr.node for tr in self.tripods[int(key)].actuatedarms]

        self.moveMotors(Key.Z,self.angleInit,self.angleMiddle,key,self.currentActuator)
        self.moveMotors(Key.E,self.angleMiddle,self.angleTop,key,self.currentActuator)
        self.moveMotors(Key.T,self.angleTop,self.angleInit,key,self.currentActuator)


        if key == Key.W :
            self.moveMotors(Key.W,self.angleInit,self.angleTop,key,[tr.node for tr in self.tripods[0].actuatedarms])
            self.moveMotors(Key.W,self.angleInit,-1.4965,key,[tr.node for tr in self.tripods[1].actuatedarms])
        if key == Key.X :
            self.moveMotors(Key.X,self.angleTop,self.angleInit,key,[tr.node for tr in self.tripods[0].actuatedarms])
            self.moveMotors(Key.X,-1.4965,self.angleInit,key,[tr.node for tr in self.tripods[1].actuatedarms])
        if key == Key.C :
            self.moveMotors(Key.C,self.angleInit,-1.4965,key,[tr.node for tr in self.tripods[0].actuatedarms])
            self.moveMotors(Key.C,self.angleInit,self.angleTop,key,[tr.node for tr in self.tripods[1].actuatedarms])
        if key == Key.V :
            self.moveMotors(Key.V,self.angleInit,self.angleTop,key,[tr.node for tr in self.tripods[0].actuatedarms])
            self.moveMotors(Key.V,self.angleInit,self.angleTop,key,[tr.node for tr in self.tripods[1].actuatedarms])



        if key == Key.uparrow:
            self.currentActuator[0].angleIn = self.currentActuator[0].angleOut + self.stepsize
        elif key == Key.downarrow:
            self.currentActuator[0].angleIn = self.currentActuator[0].angleOut - self.stepsize

        if key == Key.leftarrow:
            self.currentActuator[1].angleIn = self.currentActuator[1].angleOut + self.stepsize
        elif key == Key.rightarrow:
            self.currentActuator[1].angleIn = self.currentActuator[1].angleOut - self.stepsize

        if key == Key.plus:
            self.currentActuator[2].angleIn = self.currentActuator[2].angleOut + self.stepsize
        elif key == Key.minus:
            self.currentActuator[2].angleIn = self.currentActuator[2].angleOut - self.stepsize

class TripodController(Sofa.PythonScriptController): #Sofa.Core.Controller):
    """This controller has two roles:
       - if the user presses up/left/right/down/plus/minus, the servomotor angle
         is changed.
       - if thr user presses A, an animation is started to move the servomotor to the initial position
         of the real robot.
    """

    def __init__(self, node, tripod, actuators):
        self.stepsize = 0.01
        self.actuators = actuators
        self.node = node
        self.tripod = tripod
        print("------------------------------------")
        print(self.actuators)
        print("------------------------------------")

        # -1.2 -----> 0.87 (19degree)------> -0.54 == 38 degree actuation total
        self.angleInit = -1.2  #-1.2 #-1.4965
        self.angleMiddle = -0.87  
        self.angleTop =  -0.35 #-0.54 #-0.6 #-0.35 #-0.6 #-0.4 # -0.7 #-0.85 # -0.6

        self.dataSimu = []

    def onKeyPressed(self, key):
        self.initTripod(key)
        # self.moveMotors(Key.Z,self.angleInit,self.angleTop,key)

        if key == Key.D:
            pos = self.tripod.RigidifiedStructure.RigidParts.dofs.position
            rigid_center = pos[3]
            angle = ToEulerAngles(rigid_center[3],rigid_center[4],rigid_center[5],rigid_center[6])
            print("-----------> ANGLE CENTER SIMULATED"+self.tripod.name+": "+str(to_euler(angle[0]))+" | "+str(to_euler(angle[1]))+" | "+str(to_euler(angle[2])))


        if key == Key.Y:
            # print(spiral(self.tripod.ElasticBody.ElasticMaterialObject.container,self.tripod.ElasticBody.ElasticMaterialObject.dofs))
            print(spiral(self.tripod.ElasticBody.ElasticMaterialObject.highOrder_node.container,self.tripod.ElasticBody.ElasticMaterialObject.highOrder_node.dofs))



        self.animateTripod(key)

    def initTripod(self, key):
        if key == Key.A:
            animate(setupanimation, 
                    {"actuators": self.actuators, "step": 35.0, "angularstep": self.angleInit},
                    duration=1,
                    onDone=saveTripodPosition)

    def moveMotors(self,keyToActivate,angularstart,angularend,key,actuators,offset=None,down=False):
        if key == keyToActivate:
            animate(moveMotorSynchronously, 
                    {"actuators": actuators,"angularstart": angularstart, "angularend": angularend,'offset':offset, 'down':down}, 
                    duration=1)

    def animateTripod(self, key):

        # self.moveMotors(Key.E,self.angleTop,self.angleInit,key)

        self.moveMotors(Key.Z,self.angleInit,self.angleMiddle,key,self.actuators)
        self.moveMotors(Key.E,self.angleMiddle,self.angleTop,key,self.actuators)
        self.moveMotors(Key.T,self.angleTop,self.angleInit,key,self.actuators)


        if key == Key.uparrow:
            self.actuators[0].angleIn = self.actuators[0].angleOut + self.stepsize
        elif key == Key.downarrow:
            self.actuators[0].angleIn = self.actuators[0].angleOut - self.stepsize

        if key == Key.leftarrow:
            self.actuators[1].angleIn = self.actuators[1].angleOut + self.stepsize
        elif key == Key.rightarrow:
            self.actuators[1].angleIn = self.actuators[1].angleOut - self.stepsize

        if key == Key.plus:
            self.actuators[2].angleIn = self.actuators[2].angleOut + self.stepsize
        elif key == Key.minus:
            self.actuators[2].angleIn = self.actuators[2].angleOut - self.stepsize

######################################################################################################################
######################################################################################################################

class TripodControllerWithCom(TripodController):
    """This controller has three roles:
       - if the user presses up/left/right/down/plus/minus, the servomotor angle
         is changed.
       - if thr user presses A, an animation is started to move the servomotor to the initial position
         of the real robot.
       - if thr user presses B start the communication with controller card, send
         servomotor commands
    """

    def __init__(self, node, tripod, actuators, serialportctrl):
        TripodController.__init__(self, node, tripod, actuators)
        self.serialportctrl = serialportctrl
        # self.dofs = self.serialportctrl.dofs
        self.dofs = self.serialportctrl.Markers

        self.init = False
        self.angleCalibration = [-1.4965]*len(self.actuators)
        self.offset = [-0.3, -0.2, -0.2]

    def moveMotors(self,keyToActivate,angularstart,angularend,key,actuator,offset=None,down=False):
        if key == keyToActivate:
            # pos_real = self.dofs.position
            # rigid_center = pos_real[0]
            # angle = ToEulerAngles(rigid_center[3],rigid_center[4],rigid_center[5],rigid_center[6])

            # pos = self.dofs.position
            # rigid_center = pos[3]
            # angle = ToEulerAngles(rigid_center[3],rigid_center[4],rigid_center[5],rigid_center[6])
            # self.dataSimu.append([to_euler(angle[1]),rigid_center[1]])

            animate(moveMotorSynchronously, 
                    {"actuators": actuator,"angularstart": angularstart, "angularend": angularend,'offset':offset, 'down':down}, 
                    duration=1)

    def initTripod(self, key):
        if key == Key.A and self.serialportctrl.state == "init":
            self.serialportctrl.state = "no-comm"
            animate(setupanimation, {"actuators": self.actuators, "step": 35.0, "angularstep": self.angleInit}, duration=0.5) # -1.4965

        # Inclusion of the keystroke to start data sending = establishing communication ('comm')
        if key == Key.B and self.serialportctrl.state == "no-comm":
            self.serialportctrl.state = "comm"

    def animateTripod(self, key):
        # print("here",Key.uparrow,key)
        # self.moveMotors(Key.Z,self.angleInit,self.angleMiddle,key,self.serialportctrl,False)

        self.moveMotors(Key.Z,self.angleInit,self.angleMiddle,key,self.actuators,self.serialportctrl,True)#[self.actuators[0]],self.serialportctrl,True)
        self.moveMotors(Key.E,self.angleMiddle,self.angleTop,key,self.actuators,self.serialportctrl,True) #self.serialportctrl,True)
        self.moveMotors(Key.T,self.angleTop,self.angleInit,key, self.actuators,self.serialportctrl,True)#self.serialportctrl,True)

        # self.moveMotors(Key.E,self.angleInit,self.angleTop,key,[self.actuators[1]],self.serialportctrl,True)
        # self.moveMotors(Key.T,self.angleInit,self.angleTop,key,[self.actuators[2]],self.serialportctrl,True)
        if key == Key.Q:
            tmp =self.angleInit
            self.angleInit = self.angleTop
            self.angleTop = tmp

        if key == Key.D:
            # pos_real = self.tripod.RigidifiedStructure.optitrack.Markers.position
            # rigid_center_real = pos_real[0]
            # angle_real = ToEulerAngles(rigid_center_real[3],rigid_center_real[4],rigid_center_real[5],rigid_center_real[6])
            # print("-----------> ANGLE CENTER REAL"+self.tripod.name+": "+str(to_euler(angle_real[0]))+" | "+str(to_euler(angle_real[1]))+" | "+str(to_euler(angle_real[2])))

            print(self.serialportctrl.positions_real)
            print("################################\n################################")
            print(self.serialportctrl.positions_simu)
            print("################################\n################################")
            print(self.dataSimu)
            for i in self.dataSimu:
                print(i[1]*100)

        if key == Key.I:
            self.init = not self.init
            self.serialportctrl.calibrate = self.init 

        if not self.init:
            if key == Key.uparrow:
                self.actuators[0].angleIn = self.actuators[0].angleOut + self.stepsize
            elif key == Key.downarrow:
                self.actuators[0].angleIn = self.actuators[0].angleOut - self.stepsize

            if key == Key.leftarrow:
                self.actuators[1].angleIn = self.actuators[1].angleOut + self.stepsize
            elif key == Key.rightarrow:
                self.actuators[1].angleIn = self.actuators[1].angleOut - self.stepsize

            if key == Key.plus:
                self.actuators[2].angleIn = self.actuators[2].angleOut + self.stepsize
            elif key == Key.minus:
                self.actuators[2].angleIn = self.actuators[2].angleOut - self.stepsize
        else:

            if key == Key.uparrow:
                self.actuators[0].angleIn = self.actuators[0].angleOut + self.stepsize
            elif key == Key.downarrow:
                self.actuators[0].angleIn = self.actuators[0].angleOut - self.stepsize

            if key == Key.leftarrow:
                self.actuators[1].angleIn = self.actuators[1].angleOut + self.stepsize
            elif key == Key.rightarrow:
                self.actuators[1].angleIn = self.actuators[1].angleOut - self.stepsize

            if key == Key.plus:
                self.actuators[2].angleIn = self.actuators[2].angleOut + self.stepsize
            elif key == Key.minus:
                self.actuators[2].angleIn = self.actuators[2].angleOut - self.stepsize

            # if key == Key.uparrow:
            #     self.offset[0] += self.stepsize
            #     self.angleCalibration[0] += self.offset[0]
            # elif key == Key.downarrow:
            #     self.offset[0] -= self.stepsize
            #     self.angleCalibration[0] += self.offset[0]

            # if key == Key.leftarrow:
            #     self.offset[1] += self.stepsize
            #     self.angleCalibration[1] += self.offset[1]
            # elif key == Key.rightarrow:
            #     self.offset[1] -= self.stepsize
            #     self.angleCalibration[1] += self.offset[1]

            # if key == Key.plus:
            #     self.offset[2] += self.stepsize
            #     self.angleCalibration[2] += self.offset[2]                
            # elif key == Key.minus:
            #     self.offset[2] -= self.stepsize
            #     self.angleCalibration[2] += self.offset[2]

            # angles = []
            # for angle in self.offset:
            #     # Conversion of the angle values from radians to degrees
            #     angleDegree = angle*360/(2.0*pi)
            #     # print(angleDegree)
            #     angleByte = int(floor(angleDegree)) + 179

            #     # Limitation of the angular position's command
            #     if angleByte < 60:
            #         angleByte = 60
            #     if angleByte > 180:
            #         angleByte = 180

            #     # Filling the list of the 3 angle values
            #     angles.append(angleByte)
            # print("OFFSET :"+str(self.offset))
            # self.serialportctrl.serialport.packetOut = angles

# Data sending controller
class SerialPortController(Sofa.PythonScriptController): #Sofa.Core.Controller):
    def __init__(self, node, inputs, serialport,Markers,dofs):
        self.name = "serialportcontroller"
        self.actuatedarms = inputs
        self.serialport = serialport
        # self.serialport.packetOut = [150, 150, 150]
        self.state = "init"
        self.calibrate = False
        self.manual = [-1.4965, -1.4965, -1.4965]
        self.Markers = Markers
        self.dofs = dofs
        self.positions_real = []
        self.positions_simu = []

    def onEndAnimationStep(self, dt):
        # Data sending if the robot is initializing or in the no-communication sate
        # if self.state == "init":
        #     return

        # if self.state == "no-comm":
        #     return

        # Vector storing the simulated servomotors' angular position
        angles = []

        # self.calibrate = True

        if not self.calibrate:
            pass
            # for i,arm in enumerate(self.actuatedarms):
            #     self.manual[i] = arm.angleIn
            #     # Conversion of the angle values from radians to degrees
            #     angleDegree = arm.angleIn*360/(2.0*pi)
            #     # print(angleDegree)
            #     angleByte = int(floor(angleDegree)) + 179

            #     # Limitation of the angular position's command
            #     if angleByte < 60:
            #         angleByte = 60
            #     if angleByte > 180:
            #         angleByte = 180

            #     # Filling the list of the 3 angle values
            #     angles.append(angleByte)
            # # print("#################")
            # # The controller board of the real robot receives `angles` values
            # # print("here : "+str(angles))
            # self.serialport.packetOut = angles
        else:
            if self.Markers != None:
                pos_real = self.Markers.position
                rigid_center_real = pos_real[0]
                angle_real = ToEulerAngles(rigid_center_real[3],rigid_center_real[4],rigid_center_real[5],rigid_center_real[6])
                self.positions_real.append(angle_real)
                print("-----------> ANGLE CENTER REAL : "+str(to_euler(angle_real[0]))+" | "+str(to_euler(angle_real[1]))+" | "+str(to_euler(angle_real[2])))
            
            pos = self.dofs.position
            rigid_center = pos[3]
            angle = ToEulerAngles(rigid_center[3],rigid_center[4],rigid_center[5],rigid_center[6])
            self.positions_simu.append(angle)
            print("-----------> ANGLE CENTER SIMULATED : "+str(to_euler(angle[0]))+" | "+str(to_euler(angle[1]))+" | "+str(to_euler(angle[2])))

            for angleIn in self.manual:
                # Conversion of the angle values from radians to degrees
                angleDegree = angleIn*360/(2.0*pi)
                # print(angleDegree)
                angleByte = abs(int(floor(angleDegree))) #+ 179

                # # Limitation of the angular position's command
                # if angleByte < 60:
                #     angleByte = 60
                # if angleByte > 180:
                #     angleByte = 180

                # Filling the list of the 3 angle values
                angles.append(angleByte)

            Angles = [] #int(i) for i in angles]
            for i,arm in enumerate(self.actuatedarms):
                Angles.append(135 - abs(int(arm.angleIn*360/(2.0*pi))))

            print(Angles)   
            for i in range(3):
                if Angles[i] < 45:
                    Angles[i] = 45
                    print("----------------------->     45")
                if Angles[i] > 130: 
                    Angles[i] = 130
                    print("----------------------->     130")

                # Angles[i] = 135 - Angles[i]

            # print("#################")
            # The controller board of the real robot receives `angles` values
            print("calibrating : "+str(Angles))
            self.serialport.packetOut = Angles
