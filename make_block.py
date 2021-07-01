# make geometries for faults of varying roughness
# by changing bonds along mesh surface 

from gengeo import *
import math
import random
import sys

# x, y, z
# x= width, y =height, z= in/out of board

arg_list = sys.argv

# nohup python make_block.py blockT45L10X0Y0 1000e5 > geo5.out &
# nohup python make_block.py blockNAN 1000e5 > geo5.out &
ffile = arg_list[1] # geometry of fault included in the model, not including 'F1.txt'
nstr = arg_list[2] # number of seeds to include and distance to edge of top and bottom, 10e5, 100e5, 1000e5

print "Fault file:", ffile

is_f = 1
if 'NAN' in ffile:
    is_f =0

geofile = 'geo_'+ffile+'n'+nstr
print geofile

xd = 40.0
yd = 80.0
zd = 10.0 

# particle radii 
rmin = 0.1 ##### 0.1
rmax = 1.0

# tags
tago = 1 # off-fault
tagof = 2 # off-fault to on- fault
tagf = 3 # fault zone bonds

# maybe this should be smaller in size to allow particles to go
# from one side to the other
minPointTbl = Vector3(-5.*rmax,-5.*rmax,-5.*rmax)
maxPointTbl = Vector3(xd+5.*rmax,yd+5.*rmax,zd+5.*rmax)

# neighbour table
mntable = MNTable3D(
    minPoint=minPointTbl, 
    maxPoint=maxPointTbl,
    gridSize=2.5*rmax,
    numGroups=1)

min_pt = Vector3(0.0, 0.0, 0.0)
max_pt = Vector3(xd, yd, zd)

if is_f==1:
    print "Making Faults"
    # get the first fault
    txt = open('surfs/'+ffile+'F1.txt', 'r')

    str = txt.read()
    data = str.split("\n")
    pt3s = []
    for pt in data:
        pts = pt.split()
        if len(pt)>=3 and ('=' not in pts):
            x = float(pts[0])
            y = float(pts[1])
            z = float(pts[2])

            pt3s.append(Vector3(x, y, z))

    # make series of triangles to build into mesh
    f1pyrs = TriPatchSet()
    # make triangles out of points in file
    vi = 0
    while vi<len(pt3s):
        v1 = pt3s[vi]
        v2 = pt3s[vi+1]
        v3 = pt3s[vi+2]
        f1pyrs.addTriangle(v1, v2, v3, 1)

        vi=vi+3

    # get the second fault
    txt = open('surfs/'+ffile+'F2.txt', 'r')

    str = txt.read()
    data = str.split("\n")
    pt3s = []
    for pt in data:
        pts = pt.split()
        if len(pt)>=3 and ('=' not in pts):
            x = float(pts[0])
            y = float(pts[1])
            z = float(pts[2])

            pt3s.append(Vector3(x, y, z))

    # make series of triangles to build into mesh
    f2pyrs = TriPatchSet()
    # make triangles out of points in file
    vi = 0
    while vi<len(pt3s):
        v1 = pt3s[vi]
        v2 = pt3s[vi+1]
        v3 = pt3s[vi+2]
        f2pyrs.addTriangle(v1, v2, v3, 1)

        vi=vi+3

# iteration parameters
insertFails = 1000
maxIter = 1000
tol = 1.0e-6

#packer = InsertGenerator3D(rmin, rmax, insertFails, maxIter, tol, False) # keep same position of particles for every model
packer = InsertGenerator3D(rmin, rmax, insertFails, maxIter, tol, True) # make randomly packed particles

# make box of model
# total block volume
box = BoxWithPlanes3D(min_pt,max_pt)

# boundary planes
bottomPlane=Plane(min_pt,Vector3(0.0,1.0,0.0))
leftPlane=Plane(min_pt,Vector3(1.0,0.0,0.0))
frontPlane=Plane(min_pt,Vector3(0.0,0.0,1.0))
topPlane=Plane(max_pt,Vector3(0.0,-1.0,0.0))
rightPlane=Plane(max_pt,Vector3(-1.0,0.0,0.0))
backPlane=Plane(max_pt,Vector3(0.0,0.0,-1.0))

# add them to the box
box.addPlane(bottomPlane)
box.addPlane(leftPlane)
box.addPlane(frontPlane)
box.addPlane(topPlane)
box.addPlane(rightPlane)
box.addPlane(backPlane)

# packing full box
packer.generatePacking(volume=box,
                       ntable=mntable,
                       groupID=0,
                       tag=tago)

if is_f==1:
    # tag particles along mesh with fault tag
    mntable.tagParticlesAlongJoints(f1pyrs,
                                    0.2,
                                    tagf,
                                    -1,
                                    0)

    # tag particles along mesh with fault tag
    mntable.tagParticlesAlongJoints(f2pyrs,
                                    0.2,
                                    tagf,
                                    -1,
                                    0)

# tag random particles thoughout the model
# to try to localize deformation away from the sides
# make a sphere with 0.1 mm radius at various positions
# within y=10 20, z=0, 30, x=0, 50
# tag particles within the sphere as a fault zone
txt = open('surfs/seeds_n'+nstr+'.txt', 'r')

rad=0.5

print "adding damage spheres"
str = txt.read()
data = str.split("\n")
for pt in data:
    pts = pt.split()
    if len(pt)>=3 and ('=' not in pts):
        x = float(pts[0])
        y = float(pts[1])
        z = float(pts[2])

        print x, y, z

        #cent = Vector3(x, y, z)
        cent = Vector3(y, x, z) # switched the x and y to fit into these block models
        sphc = SphereVol(cent, rad)
        # make a text file with the location of the points
        #mntable.tagParticlesInSphere(sphc, tagf)

        mntable.tagParticlesInVolume(volume=sphc,
                                       tag = tagf,
                                       groupID = 0)


# make and tag bonds throughout off fault
mntable.generateBondsTagged(
    groupID=0,
    tolerance = 1.0e-4, 
    bondID = tago, 
    particleTag1 = tago, 
    particleTag2 = tago        
)
# make and tag bonds along fault
mntable.generateBondsTagged(
    groupID=0,
    tolerance = 1.0e-4, 
    bondID = tagf, 
    particleTag1 = tagf, 
    particleTag2 = tagf        
)
# make and tag between top and bottom block
mntable.generateBondsTagged(
    groupID=0,
    tolerance = 1.0e-4, 
    bondID = tagof, 
    particleTag1 = tago,
    particleTag2 = tagf        
)

print geofile
# write files
mntable.write(geofile+".geo",1)
mntable.write(geofile+".vtu",2)
