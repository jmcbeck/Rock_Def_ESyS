# --- geometry setup script for block with interlocking spherical particles---
from gengeo import *
from esys.lsm.geometry import *
import math
import sys

# semi-interlocking grain structure
# with walls on sides to stabilize loading

# pack volume with large spheres
# loop over spheres, increase their radius
# pack particles into overlapping spheres
# no buffer walls

# 1) make mntable for particles
# 2) make mntable for spheres (grains)

# tag of particles within cement = 0
# tag of particles in grains = 1 and above

# - input parameters --
# block dimensions
# mm
xmin = 0.0
ymin = 0.0
zmin = 0.0

# particle size range
# wider range of particle radius leads to low porosity
Rmin = 0.1 #0.1
Rmax = 1.0

arg_list = sys.argv
# range in radius of spheres before growth
Rsphere_min = float(arg_list[1])
Rsphere_max = float(arg_list[2])
growth = float(arg_list[3])
scale = float(arg_list[4])
geoItr = int(arg_list[5])

print "particle radius min and max"
print Rmin
print Rmax
print "grain radius min and max"
print Rsphere_min
print Rsphere_max
print "growth factor of grains"
print growth
print "dimension"
print scale
print "geometry iteration"
print geoItr
if geoItr>1:
        print "packing with new random seed"

filename_root = 'R'+arg_list[1]+'_'+arg_list[2]+'_'+arg_list[3]+'D'+arg_list[4]+'G'+arg_list[5]

# let user decide the scale of the model
xdim = scale
ydim = scale
zdim = scale

# ---------------------
# corner points of total domain
minPoint = Vector3(xmin,ymin,zmin)
maxPoint = Vector3(xdim,ydim,zdim)

minPointTbl = Vector3(xmin-5.*Rmax,ymin-5.*Rmax,zmin-5.*Rmax)
maxPointTbl = Vector3(xdim+5.*Rmax,ydim+5.*Rmax,zdim+5.*Rmax)

minPointTbl_sphere = Vector3(xmin-5.*Rsphere_max,ymin-5.*Rsphere_max,zmin-5.*Rsphere_max)
maxPointTbl_sphere = Vector3(xdim+5.*Rsphere_max,ydim+5.*Rsphere_max,zdim+5.*Rsphere_max)

# neighbour table for grains
mntable_spheres = MNTable3D(
        minPoint=minPointTbl_sphere, 
        maxPoint=maxPointTbl_sphere,
        gridSize=2.5*Rsphere_max,
        numGroups=1) 

# neighbour table for particles inside grains
mntable_particles = MNTable3D(
        minPoint=minPointTbl, 
        maxPoint=maxPointTbl,
        gridSize=2.5*Rmax,
        numGroups=1) # particles in grains have groupID = 0, and in walls have groupID = 1 

# -- setup packer --
# iteration parameters
insertFails = 1000
maxIter = 1000
tol = 1.0e-6

# packer
# with random seeding that changes
bool_rand = False # not random
if geoItr>1:
        bool_rand = True # random

packer_spheres = InsertGenerator3D(Rsphere_min, 
                                   Rsphere_max, 
                                   10*insertFails, 
                                   10*maxIter, 
                                   tol, 
				   bool_rand)

# box to fill with grains
minPoint = Vector3(xmin-5*Rsphere_max,ymin-5*Rsphere_max,zmin-5*Rsphere_max)
maxPoint = Vector3(xdim+5*Rsphere_max,ydim+5*Rsphere_max,zdim+5*Rsphere_max)

# total block volume
box = BoxWithPlanes3D(minPoint,maxPoint)

# boundary planes
bottomPlane=Plane(minPoint,Vector3(0.0,1.0,0.0))
leftPlane=Plane(minPoint,Vector3(1.0,0.0,0.0))
frontPlane=Plane(minPoint,Vector3(0.0,0.0,1.0))
topPlane=Plane(maxPoint,Vector3(0.0,-1.0,0.0))
rightPlane=Plane(maxPoint,Vector3(-1.0,0.0,0.0))
backPlane=Plane(maxPoint,Vector3(0.0,0.0,-1.0))

# add them to the box
box.addPlane(bottomPlane)
box.addPlane(leftPlane)
box.addPlane(frontPlane)
box.addPlane(topPlane)
box.addPlane(rightPlane)
box.addPlane(backPlane)

packer_spheres.generatePacking(volume=box, 			
                               ntable=mntable_spheres,
                               groupID=0,
                               tag=0) # tag of spheres doesn't matter 


# write a geometry file for bounding grains/spheres
print "packed the empty spheres-------------------"

packer_particles = InsertGenerator3D(Rmin, 
                                     Rmax, 
                                     insertFails, 
                                     maxIter, 
                                     tol, 
                                     False)

# pack the particles inside the small box
# box to fill with grains
minPoint = Vector3(xmin,ymin,zmin)
maxPoint = Vector3(xdim,ydim,zdim)

# total block volume
box = BoxWithPlanes3D(minPoint,maxPoint)

# boundary planes
bottomPlane=Plane(minPoint,Vector3(0.0,1.0,0.0))
leftPlane=Plane(minPoint,Vector3(1.0,0.0,0.0))
frontPlane=Plane(minPoint,Vector3(0.0,0.0,1.0))
topPlane=Plane(maxPoint,Vector3(0.0,-1.0,0.0))
rightPlane=Plane(maxPoint,Vector3(-1.0,0.0,0.0))
backPlane=Plane(maxPoint,Vector3(0.0,0.0,-1.0))

# add them to the box
box.addPlane(bottomPlane)
box.addPlane(leftPlane)
box.addPlane(frontPlane)
box.addPlane(topPlane)
box.addPlane(rightPlane)
box.addPlane(backPlane)

packer_particles.generatePacking(volume=box, 			
			ntable=mntable_particles,
			groupID=0,
			tag=0) # tag of particles not within grains

print "packed the particles------------"

# get the big grains, to pack with particles
sphere_list = mntable_spheres.getSphereListFromGroup(groupID=0)

# make list of polyhedrons in model
part_tag_i = 1 # start tag of grains in particles beyond wall particles
grain_rads = [] # keep list of the radius of all the newly created grains/spheres
grain_pos = []
for sphere in sphere_list:

        #print sphere

        center = sphere.Centre()
        rad_prev = sphere.Radius()

        xpos = center.X()
        ypos = center.Y()
        zpos = center.Z()

        # fill with particle if any part of grain within model
        if xpos > xmin-rad_prev and xpos < xdim+rad_prev and ypos > ymin-rad_prev and ypos < ydim+rad_prev and zpos > zmin-rad_prev and zpos < zdim+rad_prev:

                # grow each radius by constant factor (2)
                rad_new = rad_prev*growth
        
                sphere_pack = SphereVol(center, rad_new)

                print "taggin particles in grain--------------"
                print part_tag_i-1

                # tag particles within this sphere

                #packer_particles.generatePacking(volume=sphere_pack, 
                #                                 ntable=mntable_particles,
                #                                 groupID=part_ID,
                #                                 tag=part_tag_i)	
                mntable_particles.tagParticlesInVolume(volume=sphere_pack,
                                                       tag = part_tag_i,
                                                       groupID = 0)


                grain_rads.append(rad_new)
                grain_pos.append(center)
                part_tag_i=part_tag_i+1

max_tag = part_tag_i-1
print "Number of grains---------------"
print max_tag 

print "Removing grains with tag=0, outside grains--------------"
mntable_particles.removeParticlesWithTag(tag = 0,
                                         groupID = 0)


# bond particles in groupID=0 within grains
# particles with same tag have bondTag1 (grains, 1), particle with different have bondTag2 (between grains 0)
mntable_particles.generateClusterBonds(
        groupID= 0,
        tolerance = 1.0e-4,
        bondTag1 = 1,
        bondTag2 = 0
)


geo_filename = "geo_sandstone_"+filename_root+".geo"
vtk_filename = "vtk_sandstone_"+filename_root +".vtu"
        
print geo_filename

# write a geometry file
mntable_particles.write(geo_filename, 1)
mntable_particles.write(vtk_filename, 2)

print "Grains, radius position"
i=0
for rad in grain_rads:
    pos_curr = grain_pos[i]
    print str(pos_curr)+" "+str(rad)
    i=i+1

