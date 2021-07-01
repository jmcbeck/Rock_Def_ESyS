# triaxial_XXX.py: triaxial compression test for solid with one layer
#
#

#import the appropriate ESyS-Particle modules:
from esys.lsm import *
from esys.lsm.util import *
from esys.lsm.geometry import *
from WallLoader import WallLoaderRunnable
from ServoWallLoader import ServoWallLoaderRunnable
import math
import sys
import os

# units of length, L = 1e-3 m
# units of time, T = 1 s
# units of force, F = 1 N
# yields =>
# velocity in mm/s
# stress in MPa
# mass in 1000 kg

# set in geometry file (make_geo_ani_affin.py)
minRadius = 0.1 
maxRadius = 1.0 

#instantiate a simulation object:
sim = LsmMpi (numWorkerProcesses = 1, mpiDimList = [1,1,1])
#initialise the neighbour search algorithm:
sim.initNeighbourSearch (
    particleType = "RotSphere",
    gridSpacing = 2.5*maxRadius,
    verletDist = 2*minRadius  # MD: 0.1
)

# nohup mpirun -np 2 esysparticle run_block.py 1 50 50 1e4 1e4 1e4 geo_faultL10X0Y0 10 > runbox.out &
arg_list = sys.argv

# material parameters
coh1 = float(arg_list[1]) # cohesion of fault, 1, 10, 20, ... 50
coh2 = float(arg_list[2]) # damage zone
coh3 = float(arg_list[3]) # off-fault

E1 = float(arg_list[4]) # Young's modulus, E, of fault, 1e4=90 GPa
E2 = float(arg_list[5]) # E of damage zone
E3 = float(arg_list[6]) # E of off-fault

# geometry file
geo_str = arg_list[7] # geometry file without the .geo
confin_stress = float(arg_list[8]) # confining stress in MPa

# read in pre-made geometry file: made in make_step.py
#geo_file = geo_str+".geo"
# read in pre-made geometry file: made in make_step_nonp.py or make_step_nonp_m.py
geo_file = geo_str+".geo"

print "Using geometry file:"
print geo_file

# returns list of particles and connections
sim.readGeometry(geo_file)

# bonds between layers have same values as
# bonds within layers
fileSaveRoot = geo_str+"_C"+arg_list[1]+"_"+arg_list[2]+"_"+arg_list[3]+"_E"+arg_list[4]+"_"+arg_list[5]+"_"+arg_list[6]+"_P"+arg_list[8]

# remove files with this root, and ending with dat in directory
print "Removing pre-existing .dat files"
rm_command = 'rm *'+fileSaveRoot+'*dat'
os.system(rm_command)

print "Saving to file root:"
print fileSaveRoot

vis = 1.0
# calculate time step from max E
rho_bulk = 2700*(1e-12) # rock density
poro = 0.20 # bulk porosity of model
rho = 1e9*rho_bulk*(1-poro) 
E = max([E1, E2, E3])
k_max = E*maxRadius # max stiffness is max elasitc modulus and particle radis
mass_min = (4/3)*math.pi*rho*(minRadius**3)
dt = 0.1*math.sqrt(mass_min/k_max)
dt = round(dt, 10)
print "Time step size: ", dt
print "Density: ", rho

timeStepsTotal = 140000 

#set the number of timesteps and timestep increment:
sim.setNumTimeSteps(timeStepsTotal)
sim.setTimeStepSize(dt)

print "Cohesion:"
print "fault: ", coh1
print "damage zone: ", coh2
print "off-fault: ", coh3

print "E:"
print "fault: ", E1
print "damage zone: ", E2
print "off-fault: ", E3

# off-fault, between on- and off-fault, on-fault
bond_tags = [1, 2, 3]
cohs = [coh3, coh2, coh1]
Es = [E3, E2, E1]
intFrics = [1, 1, 1]

# dimensions of model geometry, should equal those specified in make_geo_rem_cube.py
# dimension = mm
xdim = 40.0
ydim = 80.0
zdim = 10.0

xmin = 0.0
ymin = 0.0
zmin = 0.0

# calculate force to apply to walls in order to achieve desired stress
fb_area = xdim*ydim 
rl_area = ydim*zdim 
top_area = xdim*zdim # area of top and bottom
# sigyy = fbot - ftop/area
confin_force_fb = (confin_stress*fb_area)/2.0 
confin_force_rl = (confin_stress*rl_area)/2.0 

# sigyy = (fbot - ftop)/area
# sidg = (y- (-y))/a
# sig = (2y)/a
confin_force_top = (confin_stress*top_area)/2.0 
print "Confining force on the top and bottom:", confin_force_top
print "Confining force on the left and right", confin_force_rl
print "Confining force on the front and back", confin_force_fb

# get position of top and bottom wall in geometry
minPoint = Vec3(xmin,ymin,zmin)
maxPoint = Vec3(xdim,ydim,zdim)

#create a wall at the bottom of the model:
sim.createWall (
    name = "bottom_wall",
    posn = minPoint,
    normal = Vec3(0.0000, 1.0000, 0.0000)
)

#create a wall at the top of the model:
sim.createWall (
    name = "top_wall",
    posn = maxPoint, 
    normal = Vec3(0.0000, -1.0000, 0.0000)
)

#create side wall 
sim.createWall (
    name = "left_wall",
    posn = minPoint,
    normal = Vec3(1.0000, 0.0000, 0.0000)
)

#create a wall 
sim.createWall (
    name = "right_wall",
    posn = maxPoint, 
    normal = Vec3(-1.0000, 0.0000, 0.0000)
)

#create a wall 
sim.createWall (
    name = "front_wall",
    posn = minPoint, 
    normal = Vec3(0.0000, 0.0000, 1.0000)
)

#create a wall 
sim.createWall (
    name = "back_wall",
    posn = maxPoint, 
    normal = Vec3(0.00, 0.00, -1.00)
)

# next iteration: allow different normal and shear stiffness for wall and particle interaction
# NRotSoftBondedWallPrms
# name, wallName, normalK, shearK, particleTag, tagMask, scaling (whether to scale elastic modulus by constant radius)

# add walls, then add interactions between walls and particles
#specify elastic repulsion from the bottom wall:
sim.createInteractionGroup (
    NRotElasticWallPrms (
        name = "bottom_wall_repel",
        wallName = "bottom_wall",
        normalK = E
    )
)
#specify elastic repulsion from the top wall:
sim.createInteractionGroup (
    NRotElasticWallPrms (
        name = "top_wall_repel",
        wallName = "top_wall",
        normalK = E
    )
)
#specify elastic repulsion 
sim.createInteractionGroup (
    NRotElasticWallPrms (
        name = "front_wall_repel",
        wallName = "front_wall",
        normalK = E
    )
)
#specify elastic repulsion 
sim.createInteractionGroup (
    NRotElasticWallPrms (
        name = "back_wall_repel",
        wallName = "back_wall",
        normalK = E
    )
)
#specify elastic repulsion 
sim.createInteractionGroup (
    NRotElasticWallPrms (
        name = "left_wall_repel",
        wallName = "left_wall",
        normalK = E
    )
)
#specify elastic repulsion 
sim.createInteractionGroup (
    NRotElasticWallPrms (
        name = "right_wall_repel",
        wallName = "right_wall",
        normalK = E
    )
)

#initialise frictional interactions for unbonded particles:
# possible to set friction on unbonded particles by particle tag????
# young's modulus and possion's ratio used to calculate normal and shear force on particles
sim.createInteractionGroup (
    FrictionPrms(
        name="friction",
        youngsModulus=E, # with lower E, less peak
        poissonsRatio= 0.25, 
        dynamicMu=0.1,  
        staticMu=0.1 
    )
)

# if rho not= 0, then set the particle density
if rho:
    print "Setting particle density"
    print rho
    # set particle density
    # if no mask, then all particle density should be set
    sim.setParticleDensity(
        tag = 1,
        mask = 0,
        Density = rho
    )


# loops through between layer values and within each bond
index_i = 0
for tag_i in bond_tags:
    #create rotational elastic-brittle bonds between particles:
    # bond breaks when shear stress within the bond exceeds shear strength 
    #   (cohesion+normal force*tan(internal fric angle))
    # name, youngsModulus (stress), poss, cohesion (stress), tan(angle of internal friction)
    # tag = bonds set with these values
    sim.createInteractionGroup (
        BrittleBeamPrms(
            name="bonds"+str(tag_i),
            youngsModulus= Es[index_i], 
            poissonsRatio=0.25,
            cohesion= cohs[index_i],  
            tanAngle= 1.0,
            tag=tag_i # when bond tag is this value, set the properties
        )
    )

    #create an exclusion between bonded and frictional interactions:
    sim.createExclusion (
        interactionName1 = "bonds"+str(tag_i),
        interactionName2 = "friction"
    )

    print "Setting bond tag ", tag_i
    print "Cohesion:", cohs[index_i]
    print "E:", Es[index_i]

    index_i = index_i+1


# uniaxial compression usually in quasi-static regime
# external loads applied slowle compared to compressional wavespeed
# AE generated during fracturing dissipate rapidly compared to duration of experiment
#add translational viscous damping:
# 0.002 chosen so that damping has little effect on elastic response of
# simulated rock sample, but sufficient to attenuate unwanted oscillations
sim.createInteractionGroup (
    LinDampingPrms(
        name="damping1",
        viscosity=vis, 
        maxIterations=100
    )
)
#add rotational viscous damping:
sim.createInteractionGroup (
    RotDampingPrms(
        name="damping2",
        viscosity=vis, 
        maxIterations=100
    )
)

# LOAD TOP AND BOTTOM WALLS WITH CONSTANT FORCE TO MATCH CONFINING STRESS
# ON SIDES
# only do this if confining stress > 0
if confin_stress > 0:
    print "Confining stress >  0, loading top wall with confining stress"
    # keep force on side walls constant
    # force at right points into model = neg x direction
    servo_loaderTop = ServoWallLoaderRunnable(
        LsmMpi = sim,
        interactionName = "top_wall_repel",
        force = Vec3 (0.0, -1*confin_force_top, 0.0),
        startTime = 0,
        rampTime = 40000,
        endTime = 49999 # stop constant force before displacements begin
    )
    sim.addPreTimeStepRunnable(servo_loaderTop)

    # keep force on side walls constant
    # force at front points into model = pos z direction
    servo_loaderBot = ServoWallLoaderRunnable(
        LsmMpi = sim,
        interactionName = "bottom_wall_repel", 
        force = Vec3 (0.0, confin_force_top, 0.0),
        startTime = 0,
        rampTime = 40000,
        endTime = 49999
    )
    sim.addPreTimeStepRunnable(servo_loaderBot)

# SERVO CONTROLLED SIDE WALLS WITH CONSTANT FORCE

# left side wall
servo_loaderL = ServoWallLoaderRunnable(
    LsmMpi = sim,
    interactionName = "left_wall_repel",
    force = Vec3 (confin_force_rl, 0.0, 0.0),
    startTime = 0,
    rampTime = 40000
)
sim.addPreTimeStepRunnable(servo_loaderL)

# right side wall
servo_loaderR = ServoWallLoaderRunnable(
    LsmMpi = sim,
    interactionName = "right_wall_repel", 
    force = Vec3 (-1*confin_force_rl, 0.0, 0.0),
    startTime = 0,
    rampTime = 40000
)
sim.addPreTimeStepRunnable(servo_loaderR)

# front side wall
servo_loaderF = ServoWallLoaderRunnable(
    LsmMpi = sim,
    interactionName = "front_wall_repel",
    force = Vec3 (0.0, 0.0, confin_force_fb),
    startTime = 0,
    rampTime = 40000
)
sim.addPreTimeStepRunnable(servo_loaderF)

# back side wall
servo_loaderB = ServoWallLoaderRunnable(
    LsmMpi = sim,
    interactionName = "back_wall_repel", 
    force = Vec3 (0.0, 0.0, -1*confin_force_fb),
    startTime = 0,
    rampTime = 40000
)
sim.addPreTimeStepRunnable(servo_loaderB)
 
# DISPLACE THE TOP AND BOTTOM WALL   

#add a wall loader to move the top wall:
wall_loader1 = WallLoaderRunnable(
    LsmMpi = sim,
    wallName = "top_wall",
    vPlate = Vec3 (0.0, -0.125, 0.0), 
    startTime = 50000,
    rampTime = 40000
)
sim.addPreTimeStepRunnable(wall_loader1)

# don't move bottom wall in order to match experiments
# move the bottom wall
wall_loader2 = WallLoaderRunnable(
    LsmMpi = sim,
    wallName = "bottom_wall",
    vPlate = Vec3 (0.0, 0.125, 0.0),
    startTime = 50000,
    rampTime = 40000
)
sim.addPreTimeStepRunnable(wall_loader2)

# FIELD SAVERS

#create a FieldSaver to store number of bonds:
# store number of bonds broken along layer interfaces
for tag_i in bond_tags:
    sim.createFieldSaver (
        InteractionScalarFieldSaverPrms(
            interactionName="bonds"+str(tag_i),
            fieldName="count",
            fileName="nbonds_b"+str(tag_i)+"_"+fileSaveRoot+".dat",
            fileFormat="SUM",
            beginTimeStep=0,
            endTimeStep=timeStepsTotal,
            timeStepIncr=50
        )
    )

    #create a FieldSaver to store forces between bonded particles
    #sim.createFieldSaver (
    #    InteractionVectorFieldSaverPrms(
    #        interactionName="bonds"+str(tag_i),
    #        fieldName="force",
    #        fileName="force_b"+str(tag_i)+"_"+fileSaveRoot,
    #        fileFormat="RAW_WITH_POS_ID",
    #        beginTimeStep=50000,
    #        endTimeStep=timeStepsTotal,
    #        timeStepIncr=1000
    #    )
    #)


    #create a FieldSaver to store potential energy stored in bonds:
    #sim.createFieldSaver (
    #    InteractionScalarFieldSaverPrms(
    #        interactionName="pp_bonds"+str(b_i),
    #        fieldName="potential_energy",
    #        fileName="epot_b"+str(b_i)+"_"+fileSaveRoot+".dat",
    #        fileFormat="SUM",
    #        beginTimeStep=0,
    #        endTimeStep=timeStepsTotal,
    #        timeStepIncr=50
    #    )
    #)


# NRotFriction, RotFriction
# count, sticking, slipping
# NRotBond, RotBond
# count, breaking criterion

#create a FieldSaver to store number of slipping particles, dynamic friction interactions
#sim.createFieldSaver (
#    InteractionScalarFieldSaverPrms(
#        interactionName="friction",
#        fieldName="slipping",
#        fileName="slip_"+fileSaveRoot+".dat",
#        fileFormat="RAW_WITH_POS_ID",
#        beginTimeStep=0,
#        endTimeStep=timeStepsTotal,
#        timeStepIncr=50
#    )
#)

#create a FieldSaver to store number of slipping particles, dynamic friction interactions
#sim.createFieldSaver (
#    InteractionScalarFieldSaverPrms(
#        interactionName="friction",
#        fieldName="sticking",
#        fileName="stick_"+fileSaveRoot+".dat",
#        fileFormat="RAW_WITH_POS_ID",
#        beginTimeStep=0,
#        endTimeStep=timeStepsTotal,
#        timeStepIncr=50
#    )
#)

#create a FieldSaver to store forces between unbonded particles
# sim.createFieldSaver (
#     InteractionVectorFieldSaverPrms(
#         interactionName="friction",
#         fieldName="force",
#         fileName="force_f_"+fileSaveRoot,
#         fileFormat="RAW_WITH_POS_ID",
#         beginTimeStep=50000,
#         endTimeStep=timeStepsTotal,
#         timeStepIncr=1000
#     )
# )


#create a FieldSaver to bottom and top wall positions:
posn_saver = WallVectorFieldSaverPrms(
    wallName=["bottom_wall", "top_wall"],
    fieldName="Position",
    fileName="out_position_BT_"+fileSaveRoot+".dat",
    fileFormat="RAW_SERIES",
    beginTimeStep=0,
    endTimeStep=timeStepsTotal,
    timeStepIncr=50
)
sim.createFieldSaver(posn_saver)

#create a FieldSaver for side wall positions
posn_saver_sides = WallVectorFieldSaverPrms(
    wallName=["left_wall", "right_wall", "front_wall", "back_wall"],
    fieldName="Position",
    fileName="out_position_LRFB_"+fileSaveRoot+".dat",
    fileFormat="RAW_SERIES",
    beginTimeStep=0,
    endTimeStep=timeStepsTotal,
    timeStepIncr=50
)
sim.createFieldSaver(posn_saver_sides)

#create a FieldSaver for bottom and top  wall forces:
force_saver = WallVectorFieldSaverPrms(
    wallName=["bottom_wall", "top_wall"],
    fieldName="Force",
    fileName="out_force_BT_"+fileSaveRoot+".dat",
    fileFormat="RAW_SERIES",
    beginTimeStep=0,
    endTimeStep=timeStepsTotal,
    timeStepIncr=50
)
sim.createFieldSaver(force_saver)

#create a FieldSaver to wall forces:
force_saver_sides = WallVectorFieldSaverPrms(
    wallName=["left_wall", "right_wall", "front_wall", "back_wall"], 
    fieldName="Force",
    fileName="out_force_LRFB_"+fileSaveRoot+".dat",
    fileFormat="RAW_SERIES",
    beginTimeStep=0,
    endTimeStep=timeStepsTotal,
    timeStepIncr=50
)
sim.createFieldSaver(force_saver_sides)

# unbonded growth in total kinetic energy useful if
# suspect that simulation is unstable
#create a FieldSaver to store the total kinetic energy of the particles:
#sim.createFieldSaver (
#    ParticleScalarFieldSaverPrms(
#        fieldName="e_kin",
#        fileName="ekin_"+fileSaveRoot+".dat",
#        fileFormat="SUM",
#        beginTimeStep=0,
#        endTimeStep=timeStepsTotal,
#        timeStepIncr=50
#    )
#)


# possible to check position of bonds that break as FieldSaver ?

# save checkPointer when want to see how fractures localize
# and how bonds are breaking

# create to output data
sim.createCheckPointer (
    CheckPointPrms (
        fileNamePrefix = fileSaveRoot,
        beginTimeStep = 40000,
        endTimeStep = timeStepsTotal,
        timeStepIncr = 1000
    )
)

#execute the simulation:
sim.run()
