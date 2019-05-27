# triaxial_XXX.py: triaxial compression test for solid with parallel layers
#
#

#import the appropriate ESyS-Particle modules:
from esys.lsm import *
from esys.lsm.util import *
from esys.lsm.geometry import *
from WallLoader import WallLoaderRunnable
from ServoWallLoader import ServoWallLoaderRunnable
import sys
import os

# units of length, L = 1e-3 m
# units of time, T = 1 s
# units of force, F = 1 N
# yields =>
# velocity in mm/s
# stress in MPa
# mass in 1000 kg

# set in geometry file (make_geo_rect_O.py)
minRadius = 0.1 
maxRadius = 1.0 

#instantiate a simulation object:
sim = LsmMpi (numWorkerProcesses = 1, mpiDimList = [1,1,1])
#initialise the neighbour search algorithm:
sim.initNeighbourSearch (
    particleType = "RotSphere",
    gridSpacing = 2.5*maxRadius,
    verletDist = 2*minRadius  
)

# see calc_timestep.m for calculation of time step
# time step increment = 0.1 sqrt(min particle mass/max stiffness)
# min particle mass = 4/3 pi rho Rmin^3
# max stiffness = E * Rmax 

# factor X1,000,000, Emax = 20,000, Rmin = 0.1, dt safe = 0.1
rho = 0.0018
#dt = 1.9e-6
dt = 1.0e-6 # will this help reduce wobbles?
vis = 1.0

timeStepsTotal = 500000 

#set the number of timesteps and timestep increment:
sim.setNumTimeSteps(timeStepsTotal)
sim.setTimeStepSize(dt)

arg_list = sys.argv
coh_face = float(arg_list[1])
coh_layer = float(arg_list[2])

geo_itr = float(arg_list[3])
thickness = float(arg_list[4])
theta = float(arg_list[5])
overlap = float(arg_list[6])

confin_stress = float(arg_list[7])

# read in pre-made geometry file
geo_file = 'geo_rectO_G'+arg_list[3]+'_L'+arg_list[4]+"_T"+arg_list[5]+"_O"+arg_list[6]+".geo" # Rmin = 0.1

print "Using geometry file:"
print geo_file

# returns list of particles and connections
sim.readGeometry(geo_file) # made in geo_inter_cube.py

# bonds between layers have same values as
# bonds within layers
fileSaveRoot = "rect_C"+arg_list[1]+"-"+arg_list[2]+"G"+arg_list[3]+"L"+arg_list[4]+"T"+arg_list[5]+"O"+arg_list[6]+"P"+arg_list[7]

# remove files with this root, and ending with dat in directory
#print "Removing pre-existing .dat files"
#rm_command = 'rm *'+fileSaveRoot+'*dat'
#os.system(rm_command)

print "Saving to file root:"
print fileSaveRoot

# remove existing .dat files from directory???
# otherwise will end up writing to end?

# constant bond parameters, same within layers and between layers
E_interlay = 20000.0 #50000.0 #15000.0 # higher E, yeilds more pronounced peak? wiht 10,000 E= 41, highest E = 60, cluster near 40
# micro elastic modulus within layers
E_lay1 = E_interlay 
# E of friction interactions
E_fric = E_interlay
# stiffness of wall
E_wall = E_lay1

# internal friction coeff between layers
muoI_interlay = 1.0 
# internal friction coefficient within layers
muoI_lay1 = 1.0 

# cohesion within and between layers
# changing
coh_lay1 = coh_layer 
coh_interlay = coh_face  

print "Cohesion within layers, and between layers:"
print coh_lay1
print coh_interlay

# first value = property between (inter) layers
# remaining values = property within each layer
cohs = [coh_interlay, coh_lay1, coh_lay1, coh_lay1, coh_lay1, coh_lay1, coh_lay1, coh_lay1, coh_lay1, coh_lay1, coh_lay1, coh_lay1, coh_lay1, coh_lay1, coh_lay1, coh_lay1, coh_lay1, coh_lay1, coh_lay1, coh_lay1]
Es = [E_interlay, E_lay1, E_lay1, E_lay1, E_lay1, E_lay1, E_lay1, E_lay1, E_lay1, E_lay1, E_lay1, E_lay1, E_lay1, E_lay1, E_lay1, E_lay1, E_lay1, E_lay1, E_lay1, E_lay1, E_lay1, E_lay1, E_lay1, E_lay1]
intFrics = [muoI_interlay, muoI_lay1, muoI_lay1, muoI_lay1, muoI_lay1, muoI_lay1, muoI_lay1, muoI_lay1,  muoI_lay1, muoI_lay1, muoI_lay1, muoI_lay1, muoI_lay1, muoI_lay1, muoI_lay1,  muoI_lay1, muoI_lay1, muoI_lay1, muoI_lay1, muoI_lay1,  muoI_lay1, muoI_lay1]

# bondTag = 0 when between layers
# bondTag > 0 when inside layer

# eventually read these in from geometry file
# dimension = mm
xdim = 20.0
ydim = 40.0
zdim = 20.0

xmin = 0.0
ymin = 0.0
zmin = 0.0

# calculate force to apply to walls in order to achieve desired stress
sidewall_area = xdim*ydim # area of sidewalls when xdim = zdim
top_area = xdim*zdim # area of top and bottom
confin_force_side = 0.5*confin_stress*sidewall_area 
confin_force_top = 0.5*confin_stress*top_area 

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
        normalK = E_wall
    )
)
#specify elastic repulsion from the top wall:
sim.createInteractionGroup (
    NRotElasticWallPrms (
        name = "top_wall_repel",
        wallName = "top_wall",
        normalK = E_wall
    )
)
#specify elastic repulsion 
sim.createInteractionGroup (
    NRotElasticWallPrms (
        name = "front_wall_repel",
        wallName = "front_wall",
        normalK = E_wall
    )
)
#specify elastic repulsion 
sim.createInteractionGroup (
    NRotElasticWallPrms (
        name = "back_wall_repel",
        wallName = "back_wall",
        normalK = E_wall
    )
)
#specify elastic repulsion 
sim.createInteractionGroup (
    NRotElasticWallPrms (
        name = "left_wall_repel",
        wallName = "left_wall",
        normalK = E_wall
    )
)
#specify elastic repulsion 
sim.createInteractionGroup (
    NRotElasticWallPrms (
        name = "right_wall_repel",
        wallName = "right_wall",
        normalK = E_wall
    )
)

#initialise frictional interactions for unbonded particles:
# possible to set friction on unbonded particles by particle tag????
# young's modulus and possion's ratio used to calculate normal and shear force on particles
sim.createInteractionGroup (
    FrictionPrms(
        name="friction",
        youngsModulus=E_fric, # with lower E, less peak
        poissonsRatio= 0.25, 
        dynamicMu=0.1,  
        staticMu=0.1 
    )
)

particleList = sim.getParticleList()
part_tags = []
for part in particleList:
    part_tag = part.getTag()
    part_tags.append(part_tag)

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


part_max = int(max(part_tags))
bond_max = int(max(part_tags))+1

print "Max number of bond tags"
print bond_max

# loops through between layer values and within each bond
tag_i = 0
while tag_i <= bond_max:
    #create rotational elastic-brittle bonds between particles:
    # bond breaks when shear stress within the bond exceeds shear strength 
    #   (cohesion+normal force*tan(internal fric angle))
    # name, youngsModulus (stress), poss, cohesion (stress), tan(angle of internal friction)
    # tag = bonds set with these values
    sim.createInteractionGroup (
        BrittleBeamPrms(
            name="pp_bonds"+str(tag_i),
            youngsModulus= Es[tag_i], 
            poissonsRatio=0.25,
            cohesion= cohs[tag_i],  
            tanAngle= intFrics[tag_i], 
            tag=tag_i # when bond tag is this value, set the properties
        )
    )

    #create an exclusion between bonded and frictional interactions:
    sim.createExclusion (
        interactionName1 = "pp_bonds"+str(tag_i),
        interactionName2 = "friction"
    )

    tag_i = tag_i+1


# uniaxial compression usually in quasi-static regime
# external loads applied slowle compared to compressional wavespeed
# AE generated during fracturing dissipate rapidly compared to duration of experiment
#add translational viscous damping:
# 0.002 chosen so that damping has little effect on elastic response of
# simulated rock sample, but sufficient to attenuate unwanted oscillations
sim.createInteractionGroup (
    LinDampingPrms(
        name="damping1",
        viscosity=vis, # for quasi static as large as 0.5 0.1,  0.002
        maxIterations=100
    )
)
#add rotational viscous damping:
sim.createInteractionGroup (
    RotDampingPrms(
        name="damping2",
        viscosity=vis, # 0.1
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
        interactionName = "top_wall_repel", #"rwall_bonds",
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
        interactionName = "bottom_wall_repel", #"fwall_bonds",
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
    force = Vec3 (confin_force_side, 0.0, 0.0),
    startTime = 0,
    rampTime = 40000
)
sim.addPreTimeStepRunnable(servo_loaderL)

# right side wall
servo_loaderR = ServoWallLoaderRunnable(
    LsmMpi = sim,
    interactionName = "right_wall_repel", 
    force = Vec3 (-1*confin_force_side, 0.0, 0.0),
    startTime = 0,
    rampTime = 40000
)
sim.addPreTimeStepRunnable(servo_loaderR)

# front side wall
servo_loaderF = ServoWallLoaderRunnable(
    LsmMpi = sim,
    interactionName = "front_wall_repel",
    force = Vec3 (0.0, 0.0, confin_force_side),
    startTime = 0,
    rampTime = 40000
)
sim.addPreTimeStepRunnable(servo_loaderF)

# back side wall
servo_loaderB = ServoWallLoaderRunnable(
    LsmMpi = sim,
    interactionName = "back_wall_repel", 
    force = Vec3 (0.0, 0.0, -1*confin_force_side),
    startTime = 0,
    rampTime = 40000
)
sim.addPreTimeStepRunnable(servo_loaderB)
 
# DISPLACE THE TOP AND BOTTOM WALL   

#add a wall loader to move the top wall:
wall_loader1 = WallLoaderRunnable(
    LsmMpi = sim,
    wallName = "top_wall",
    vPlate = Vec3 (0.0, -0.125, 0.0), # -0.125
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
b_i=0
while b_i <= bond_max:
    sim.createFieldSaver (
        InteractionScalarFieldSaverPrms(
            interactionName="pp_bonds"+str(b_i),
            fieldName="count",
            fileName="nbonds_b"+str(b_i)+"_"+fileSaveRoot+".dat",
            fileFormat="SUM",
            beginTimeStep=0,
            endTimeStep=timeStepsTotal,
            timeStepIncr=50
        )
    )

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

    b_i=b_i+1

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
#sim.createCheckPointer (
#    CheckPointPrms (
#        fileNamePrefix = fileSaveRoot,
#        beginTimeStep = 0,
#        endTimeStep = timeStepsTotal,
#        timeStepIncr = 100
#    )
#)

#execute the simulation:
sim.run()
