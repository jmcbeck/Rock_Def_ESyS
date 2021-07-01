# make geometries for faults of varying roughness
# by changing bonds along mesh surface 

from gengeo import *
import math
import random
import sys

# x, y, z
# x= long dim, y =height, z= in/out of board

arg_list = sys.argv

# nohup python make_step_nonp_mesh.py faultL10X0Y0 1000e5 mesh5 > geo5.out &
ffile = arg_list[1] # geometry of fault included in the model, not including 'F1.txt'
nstr = arg_list[2] # number of seeds to include and distance to edge of top and bottom, 100e5, 100e10, 1000e5, 1000e10
meshr = arg_list[3] # mesh file root, mesh15, mesh10, mesh5

print "Fault file:", ffile
print "mesh root", meshr

mstr = meshr.replace('mesh', '')
#geofile = 'geo_'+ffile+'n'+nstr+'m'+mstr
# added B for big
geofile = 'geo_B'+ffile+'n'+nstr+'m'+mstr
print geofile

xd = 80.0
yd = 40.0
zd = 10.0 

# assume if we changed it 80 mm tall
if 'DY' in mstr:
    yd = 80.0

# particle radii 
rmin = 0.1 ##### 0.1
rmax = 1.0

# tags
tago = 1 # off-fault
tagof = 2 # off-fault to on- fault
tagf = 3 # fault zone bonds
tagtop = 4 # tag of particles along the top
tagbot = 5 # tag of particles along the bottom
taglt = 6 # tag of particle on left top
taglb = 7 # on left bottom
tagrt = 8 # tag of particle on left top
tagrb = 9 # on left bottom

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

# neighbour table
# for periodic boundaries in 3D
# mntable = CircMNTable3D(
#     minPoint=min_pt, 
#     maxPoint=max_pt,
#     gridSize=2.5*rmax,
#     numGroups=1)



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


# make the mesh walls on the top right
txt = open('surfs/'+meshr+'_rt.msh', 'r')
str = txt.read()
data = str.split("\n")
pt3s = []
wall_rt = TriPatchSet()
for pt in data:
    pts = pt.split()

    if len(pts)==6:
        x = float(pts[3])
        y = float(pts[4])
        z = float(pts[5])
        pt3s.append(Vector3(x, y, z))

    if len(pts)==5:
        i1 = int(pts[2])
        i2 = int(pts[3])
        i3 = int(pts[4])

        v1 = pt3s[i1]
        v2 = pt3s[i2]
        v3 = pt3s[i3]
        wall_rt.addTriangle(v1, v2, v3, 1)      


# make the mesh walls on the bottom right
txt = open('surfs/'+meshr+'_rb.msh', 'r')
str = txt.read()
data = str.split("\n")
pt3s = []
wall_rb = TriPatchSet()
for pt in data:
    pts = pt.split()

    if len(pts)==6:
        x = float(pts[3])
        y = float(pts[4])
        z = float(pts[5])
        pt3s.append(Vector3(x, y, z))

    if len(pts)==5:
        i1 = int(pts[2])
        i2 = int(pts[3])
        i3 = int(pts[4])

        v1 = pt3s[i1]
        v2 = pt3s[i2]
        v3 = pt3s[i3]
        wall_rb.addTriangle(v1, v2, v3, 1)  

# make the mesh walls on the top left
txt = open('surfs/'+meshr+'_lt.msh', 'r')
str = txt.read()
data = str.split("\n")
pt3s = []
wall_lt = TriPatchSet()
for pt in data:
    pts = pt.split()

    if len(pts)==6:
        x = float(pts[3])
        y = float(pts[4])
        z = float(pts[5])
        pt3s.append(Vector3(x, y, z))

    if len(pts)==5:
        i1 = int(pts[2])
        i2 = int(pts[3])
        i3 = int(pts[4])

        v1 = pt3s[i1]
        v2 = pt3s[i2]
        v3 = pt3s[i3]
        wall_lt.addTriangle(v1, v2, v3, 1) 

# make the mesh walls on the bottom left
txt = open('surfs/'+meshr+'_lb.msh', 'r')
str = txt.read()
data = str.split("\n")
pt3s = []
wall_lb = TriPatchSet()
for pt in data:
    pts = pt.split()

    if len(pts)==6:
        x = float(pts[3])
        y = float(pts[4])
        z = float(pts[5])
        pt3s.append(Vector3(x, y, z))

    if len(pts)==5:
        i1 = int(pts[2])
        i2 = int(pts[3])
        i3 = int(pts[4])

        v1 = pt3s[i1]
        v2 = pt3s[i2]
        v3 = pt3s[i3]
        wall_lb.addTriangle(v1, v2, v3, 1) 

# iteration parameters
insertFails = 1000
maxIter = 1000
tol = 1.0e-6

#packer = InsertGenerator3D(rmin, rmax, insertFails, maxIter, tol, False) # keep same position of particles for every model
packer = InsertGenerator3D(rmin, rmax, insertFails, maxIter, tol, True) # make randomly packed particles
# add indication of periodic boundaries here?
# circDimList = [True, False, False],

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
        #cents.append(Vector3(x, y, z))
        cent = Vector3(x, y, z)
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

# tag the particles on the top right, top left, bottom right and bottom left
# tag particles along mesh with fault tag
mntable.tagParticlesAlongJoints(wall_lt, 1.0, taglt, -1, 0)
mntable.tagParticlesAlongJoints(wall_lb, 1.0, taglb, -1, 0)
mntable.tagParticlesAlongJoints(wall_rt, 1.0, tagrt, -1, 0)
mntable.tagParticlesAlongJoints(wall_rb, 1.0, tagrb, -1, 0)

# tag the top and bottom planes with particles in order to move them
mntable.tagParticlesAlongPlane(topPlane, 1.0, tagtop)
mntable.tagParticlesAlongPlane(bottomPlane, 1.0, tagbot)

print geofile
# write files
mntable.write(geofile+".geo",1)
mntable.write(geofile+".vtu",2)
