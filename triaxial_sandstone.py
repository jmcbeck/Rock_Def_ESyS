# triaxial_sandstone.py: triaxial compression test for sandstones
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

# set in geometry file
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

# porosity=30%, density factor X1e6, Emax = 20,000, dt safe = 0.1
# rmin, rmax = 0.1, 1.0
rho = 0.0019
dt = 1.2e-6
vis = 0.001 # 0.002 in tutorial, influencing failure...YES!

timeStepsTotal = 300000 # possible to reduce?
#timeStepsTotal = 1000000 

#set the number of timesteps and timestep increment:
sim.setNumTimeSteps(timeStepsTotal)
sim.setTimeStepSize(dt)

arg_list = sys.argv
coh_inter = float(arg_list[1])
coh_grains = float(arg_list[2])

grain_Rmin = float(arg_list[3])
grain_Rmax = float(arg_list[4])
grain_growth = float(arg_list[5])
dim = float(arg_list[6])

geoI = float(arg_list[7])
confin_stress = float(arg_list[8])

# read in pre-made geometry file
geo_file = 'geo_sandstone'+'_R'+arg_list[3]+'_'+arg_list[4]+"_"+arg_list[5]+"D"+arg_list[6]+"G"+arg_list[7]+".geo"

# returns list of particles and connections
sim.readGeometry(geo_file) # made in make_geo_grains.py
print "Using geometry file:"
print geo_file

# bonds between layers have same values as
# bonds within layers
fileSaveRoot = "triSANDSTONE_C"+arg_list[1]+"_"+arg_list[2]+"R"+arg_list[3]+"_"+arg_list[4]+"_"+arg_list[5]+"D"+arg_list[6]+"G"+arg_list[7]+"P"+arg_list[8]

print "Saving to file root:"
print fileSaveRoot

# constant bond elastic modulus parameters
E_inter = 50000.0 # 100,000.0 in tutorial, layers = 20,000.0
E_grains = E_inter
# E of friction interactions
E_fric = E_inter
# stiffness of moving walls
E_wall = E_inter

# bonds between multiple grains (0), bonds within grains (1)
bond_tags = [0, 1] 

# bonds IDs = 0) grain-grain, grain-cement and w/in cement, (1) w/in grains
# should be low, high, low
cohs = [coh_inter,  coh_grains]
Es = [E_inter, E_grains]

# eventually read these in from geometry file
# dimension = mm
xdim= 20.0 #dim does this matter??
ydim= 20.0 #dim
zdim= 20.0 #dim

xmin = 0.0
ymin = 0.0
zmin = 0.0

sidewall_area = xdim*ydim # area of sidewalls when xdim = zdim
top_area = xdim*zdim # area of top and bottom
confin_force = 0.5*confin_stress*sidewall_area # in N, applied on sides
confin_force_top = 0.5*confin_stress*top_area # axial stress calculated from equal oppo directions, make 1/2 to match sides

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

#specify elastic repulsion 
sim.createInteractionGroup (
    NRotElasticWallPrms (
        name = "bottom_wall_repel",
        wallName = "bottom_wall",
        normalK = E_wall
    )
)
#specify elastic repulsion 
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
        youngsModulus=E_fric, 
        poissonsRatio= 0.25, 
        dynamicMu=0.001, #0.4, in tutorial
        staticMu=0.001 #0.6
    )
)

# if rho not= 0, then set the particle density
if rho:
    print "Setting particle density"
    print rho
    # set particle density
    # if no mask, then all particle density should be set ??
    sim.setParticleDensity(
        tag = 1,
        mask = 0,
        Density = rho
    )


# loops through list of bond numbers in model [0, 1, 10, 100]
index_i=0
for tag_i in bond_tags:

    #create rotational elastic-brittle bonds between particles:
    # bond breaks when shear stress within the bond exceeds shear strength 
    #   (cohesion+normal force*tan(internal fric angle))
    # name, youngsModulus (stress), poss, cohesion (stress), tan(angle of internal friction)
    # tag = bonds set with these values
    sim.createInteractionGroup (
        BrittleBeamPrms(
            name="pp_bonds"+str(tag_i),
            youngsModulus= Es[index_i], 
            poissonsRatio=0.25,
            cohesion= cohs[index_i],  
            tanAngle= 1.0, # 1.0
            tag=tag_i # when bond tag is this value, set the properties
        )
    )

    print "Bond tag"
    print tag_i
    print "Setting cohesion to" 
    print cohs[index_i]
    print "Setting E to"
    print Es[index_i]

    #create an exclusion between bonded and frictional interactions:
    sim.createExclusion (
        interactionName1 = "pp_bonds"+str(tag_i),
        interactionName2 = "friction"
    )

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

# only do this if confining stress > 0
# add force control initially to top and bottom wall
if confin_stress > 0:
    print "Confining stress>  0"
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
        interactionName = "bottom_wall_repel", 
        force = Vec3 (0.0, confin_force_top, 0.0),
        startTime = 0,
        rampTime = 40000,
        endTime = 49999
    )
    sim.addPreTimeStepRunnable(servo_loaderBot)

    # keep force on side walls constant
    # force at left points into model = positive
    servo_loaderL = ServoWallLoaderRunnable(
        LsmMpi = sim,
        interactionName = "left_wall_repel",
        force = Vec3 (confin_force, 0.0, 0.0),
        startTime = 0,
        rampTime = 40000
    )
    sim.addPreTimeStepRunnable(servo_loaderL)

    # keep force on side walls constant
    # force at right points into model = neg x direction
    servo_loaderR = ServoWallLoaderRunnable(
        LsmMpi = sim,
        interactionName = "right_wall_repel",
        force = Vec3 (-1*confin_force, 0.0, 0.0),
        startTime = 0,
        rampTime = 40000
    )
    sim.addPreTimeStepRunnable(servo_loaderR)

    # keep force on side walls constant
    # force at front points into model = pos z direction
    servo_loaderF = ServoWallLoaderRunnable(
        LsmMpi = sim,
        interactionName = "front_wall_repel", 
        force = Vec3 (0.0, 0.0, confin_force),
        startTime = 0,
        rampTime = 40000
    )
    sim.addPreTimeStepRunnable(servo_loaderF)
    
    # keep force on side walls constant
    # force at back points into model = neg z direction
    servo_loaderB = ServoWallLoaderRunnable(
        LsmMpi = sim,
        interactionName = "back_wall_repel", 
        force = Vec3 (0.0, 0.0, -1*confin_force),
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
for tag_i in bond_tags:
    sim.createFieldSaver (
        InteractionScalarFieldSaverPrms(
            interactionName="pp_bonds"+str(tag_i),
            fieldName="count",
            fileName="nbonds_b"+str(tag_i)+"_"+fileSaveRoot+".dat",
            fileFormat="SUM",
            beginTimeStep=0,
            endTimeStep=timeStepsTotal,
            timeStepIncr=100
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
    #        timeStepIncr=10
    #    )
    #)

    #create a FieldSaver to store forces between particles
    #sim.createFieldSaver (
    #    InteractionVectorFieldSaverPrms(
    #        interactionName="pp_bonds"+str(b_i),
    #        fieldName="normal_force",
    #        fileName="forceN_b"+str(b_i)+"_"+fileSaveRoot+".dat",
    #        fileFormat="RAW_WITH_POS_ID",
    #        beginTimeStep=0,
    #        endTimeStep=timeStepsTotal,
    #        timeStepIncr=10
    #    )
    #)

    #create a FieldSaver to store forces between particles
    #sim.createFieldSaver (
    #    InteractionVectorFieldSaverPrms(
    #        interactionName="pp_bonds"+str(b_i),
    #        fieldName="force",
    #        fileName="force_b"+str(b_i)+"_"+fileSaveRoot+".dat",
    #        fileFormat="RAW_WITH_POS_ID",
    #        beginTimeStep=0,
    #        endTimeStep=timeStepsTotal,
    #        timeStepIncr=10
    #    )
    #)



#create a FieldSaver to store forces between particles
#sim.createFieldSaver (
#    InteractionVectorFieldSaverPrms(
#        interactionName="friction",
#        fieldName="force",
#        fileName="force_friction_"+fileSaveRoot+".dat",
#        fileFormat="RAW_WITH_POS_ID",
#        beginTimeStep=0,
#        endTimeStep=timeStepsTotal,
#        timeStepIncr=10
#    )
#)

#create a FieldSaver to store number of slipping particles, dynamic friction interactions
#sim.createFieldSaver (
#    InteractionScalarFieldSaverPrms(
#        interactionName="friction",
#        fieldName="slipping",
#        fileName="slip_friction_"+fileSaveRoot+".dat",
#        fileFormat="RAW_WITH_POS_ID",
#        beginTimeStep=0,
#        endTimeStep=timeStepsTotal,
#        timeStepIncr=10
#    )
#)

#create a FieldSaver to store number of slipping particles, dynamic friction interactions
#sim.createFieldSaver (
#    InteractionScalarFieldSaverPrms(
#        interactionName="friction",
#        fieldName="sticking",
#        fileName="stick_friction_"+fileSaveRoot+".dat",
#        fileFormat="RAW_WITH_POS_ID",
#        beginTimeStep=0,
#        endTimeStep=timeStepsTotal,
#        timeStepIncr=10
#    )
#)


# NRotFriction, RotFriction
# count, sticking, slipping
# NRotBond, RotBond
# count, breaking criterion


#create a FieldSaver to bottom and top wall positions:
posn_saver = WallVectorFieldSaverPrms(
    wallName=["bottom_wall", "top_wall"],
    fieldName="Position",
    fileName="out_position_BT_"+fileSaveRoot+".dat",
    fileFormat="RAW_SERIES",
    beginTimeStep=0,
    endTimeStep=timeStepsTotal,
    timeStepIncr=100
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
    timeStepIncr=100
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
    timeStepIncr=100
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
    timeStepIncr=100
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
#        timeStepIncr=10
#    )
#)

# save checkPointer when want to see how fractures localize
# and how bonds are breaking

# create to output data
#sim.createCheckPointer (
#    CheckPointPrms (
#        fileNamePrefix = fileSaveRoot,
#        beginTimeStep = 5000,
#        endTimeStep = timeStepsTotal,
#        timeStepIncr = 5000
#    )
#)

#execute the simulation:
sim.run()

